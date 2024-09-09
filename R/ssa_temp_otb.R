# OTB Temperature signal processing with SSA
#
rm(list=ls(all=TRUE)) 

# Load libraries
if(!require(Rssa)) { install.packages('Rssa') }; library(Rssa)
if(!require(lubridate)) { install.packages('lubridate') }; library(lubridate)
if(!require(dplyr)) { install.packages('dplyr') }; library(dplyr)

# Load & process data
load( "../data-clean/epcwq_clean.Rdata")
# Aggregate to monthly timeframe
tempdat1 <- epcwq3[ which(epcwq3$param=="Temp_top"), c('date','value') ]
tempdat1$month <- floor_date( tempdat1$date, 'month' )
tempdat2 <- tempdat1 |> group_by(month) |>
              summarise( value = mean(value) ) |> as.data.frame()
# Create time series
tempdat <- data.frame( month = seq.Date( min(tempdat1$month),
                                         max(tempdat1$month), 'month' ) )
tempdat <- left_join( tempdat, tempdat2, by = 'month' )
# Impute missing values (linear interpolation)
which(is.na(tempdat$value))
tempdat$value[244] <- tempdat$value[243] + diff(tempdat$value[c(243,246)])/3
tempdat$value[245] <- tempdat$value[246] - diff(tempdat$value[c(243,246)])/3

# Plot time series
plot( value ~ month, data = tempdat, las = 1,
      type = 'l', lwd = 2, col = rgb(0,0,0,0.6), ylim = c(0,max(tempdat$value)),
      main = "Temp ts", xlab = '', ylab = 'Temperature (C)' )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
abline( v = seq.Date( min(tempdat$month),
                      max(tempdat$month)+60, 'year' ),
        col = rgb(0,0,0,0.1) )
points( value ~ month, data = tempdat, pch = 21, lwd = 2, cex = 0.5,
        col = rgb(0,0,0,0.6), bg = rgb(1,1,1,1) )

# Plot frequencies
# select column
colnames( tempdat )
var <- 'value'
x <- tempdat[,var] |> na.omit()
x <- x - mean(x)
# Fourier transform
spec <- spectrum( x, method = 'pgram', plot = FALSE )
df <- data.frame( power = spec$spec, period = 1/spec$freq )
df <- df[ order( df$period ), ]
# plot time series
par(mfrow=c(2,1))
plot( x, type = 'l', col = rgb(0,0,0,0.5), lwd = 3,
      las = 1, xlab = 'time', ylab = 'x',
      main = paste0( "Time series (",var,")" ),
      cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3 )
points( x, col = rgb(0,0,0,1), pch = 16, cex = 0.8 )
abline( v = axTicks(1), col = rgb(0,0,0,0.1))
abline( h = axTicks(2), col = rgb(0,0,0,0.1))  
# plot periodogram
plot( power ~ period, data = df,
      type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
      main = 'Periodogram', las = 1, bty = "L",
      cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3 )
points( power ~ period, data = df, pch = 16, cex = 0.8 )
abline( v = axTicks(1), col = rgb(0,0,0,0.1))
abline( h = axTicks(2), col = rgb(0,0,0,0.1))
df <- df[ order( df$power, decreasing = TRUE ), ]
text( x = df$period[1:2], y = df$power[1:2],
      labels = round(df$period[1:2],2), pos = 4, font = 2 )
# print spectrum
df |> head(10)


# SSA decomposition
##
# Set window length
length( x )  
win <- 144
# Decompose
obj <- ssa( x, L = win, neig = win,
            kind = 'toeplitz-ssa' )
# Eigentriple plots
# Singular values
par(mfrow=c(1,1))
obj$sigma |> plot( main = "Singular values",
                   xlab = 'eigentriple', ylab = 'singular value' )
obj$sigma |> lines()
abline( v = axTicks(1), col = rgb(0,0,0,0.1))
abline( h = axTicks(2), col = rgb(0,0,0,0.1))
# Eigenvectors
obj |> plot( type = 'vectors', numvectors = 20 )
obj |> plot( type = 'paired', numvectors = 20 )
# W-correlation matrix
obj |> plot( type = 'wcor' )
wcor(obj,groups = 1:30) |> plot()


# SSA grouping
##
# Specify signal component groups
grp <- list( c(1,2), c(3,4), c(5,8) )

# SSA reconstruction
##
# Add a residuals (noise) group
grp[[ length(grp)+1 ]] <- which( !(1:win %in% unlist(grp)) )
# Reconstruct grouped components
recon <- obj |> reconstruct( groups = grp )
# W-correlation matrix
wcor.recon <- wcor( obj, groups = grp )
wcor.recon |> plot()
wcor.recon
# Compute variance explained by each component
eigenvals <- obj$sigma^2
varexp <- lapply( grp, function(x) sum( eigenvals[x] ) * 100 / sum(eigenvals) )

# Plot reconstructed components
par(mfrow=c(3,2))
ylims <- range(unlist(recon))
for( i in 1:length(recon) ){
  if( i < length(grp) ){
    # Plot signal components
    plot( recon[[i]], type = 'l', ylim = ylims, las = 1,
          main = paste0("Group ",paste(grp[[i]],collapse=", "),
                        " (",round(varexp[[i]],2),"%)"),
          xlab = '', ylab = ''
    )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
  } else {
    # Plot noise
    plot( recon[[i]], type = 'l', col = rgb(0.1,0.2,1,0.8),
          ylim = ylims, las = 1,
          main = paste0("Noise"," (",round(varexp[[i]],2),"%)"),
          xlab = '', ylab = ''
    )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
  }
}  # // end i 

# Reconstruct signal and noise
signal <- do.call( cbind, recon[ 1:(length(recon)-1) ] ) |> rowSums()
sigstrength <- varexp[1:(length(varexp)-1)] |> unlist() |> sum()
noise <- recon[[ length(recon) ]]

# Plot signal with noise
par(mfrow=c(2,1))
ylims2 <- range( x, signal, noise )
plot( x, type = 'l', ylim = ylims2, las = 1,
      main = paste0("Signal (",round(sigstrength,2),"%)"),
      xlab = "",
      lwd = 2, col = rgb(0,0,0,0.6) )
lines( signal, lwd = 4, col = rgb(1,0.2,0.1,0.8) )
abline( v = axTicks(1), col = rgb(0,0,0,0.1))
abline( h = axTicks(2), col = rgb(0,0,0,0.1))
plot( noise, type = 'l', ylim = ylims2, las = 1,
      main = "Noise", xlab = "time",
      col = rgb(0.1,0.2,1,0.8) )
abline( v = axTicks(1), col = rgb(0,0,0,0.1))
abline( h = axTicks(2), col = rgb(0,0,0,0.1))


# Export ts, signal, noise
ssa_temp_otb <- list( dat = data.frame( month = tempdat$month,
                                       ts = tempdat$value,
                                       signal = signal + mean(tempdat$value),
                                       noise = noise + mean(tempdat$value) ),
                     grp = grp,
                     sigstrength = sigstrength 
)
save( ssa_temp_otb, file = "../data/ssa_temp_otb.RData" )
