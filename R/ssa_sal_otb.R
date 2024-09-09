# OTB top salinity signal processing with SSA
#
rm(list=ls(all=TRUE)) 

# Load libraries
if(!require(Rssa)) { install.packages('Rssa') }; library(Rssa)
if(!require(lubridate)) { install.packages('lubridate') }; library(lubridate)
if(!require(dplyr)) { install.packages('dplyr') }; library(dplyr)

# Load & process data
load( "../data-clean/epcwq_clean.Rdata")
# Aggregate to monthly timeframe
saldat1 <- epcwq3[ which(epcwq3$param=="Sal_top"), c('date','value') ]
saldat1$month <- floor_date( saldat1$date, 'month' )
saldat2 <- saldat1 |> group_by(month) |>
             summarise( value = mean(value) ) |> as.data.frame()
# Create time series
saldat <- data.frame( month = seq.Date( min(saldat1$month),
                                        max(saldat1$month), 'month' ) )
saldat <- left_join( saldat, saldat2, by = 'month' )
# Impute missing values (linear interpolation)
missing.idx <- which(is.na(saldat$value))
plot( value ~ month, data = saldat )
abline( v = saldat$month[missing.idx], col = rgb(1,0,0,0.5) )
saldat$value[78] <- mean( saldat$value[c(77,79)] )
saldat$value[86] <- mean( saldat$value[c(85,87)] )
saldat$value[114] <- saldat$value[113] + diff(saldat$value[c(113,116)])/3
saldat$value[115] <- saldat$value[116] - diff(saldat$value[c(113,116)])/3
saldat$value[122] <- mean( saldat$value[c(121,123)] )
saldat$value[133] <- mean( saldat$value[c(132,134)] )
saldat$value[244] <- saldat$value[243] + diff(saldat$value[c(243,246)])/3
saldat$value[245] <- saldat$value[246] - diff(saldat$value[c(243,246)])/3
points( value ~ month, data = saldat[missing.idx,], col = 2 )

# Plot time series
plot( value ~ month, data = saldat, las = 1,
      type = 'l', lwd = 2, col = rgb(0,0,0,0.6),
      main = "Top salinity ts", xlab = '', ylab = 'Salinity (PSU)' )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
abline( v = seq.Date( min(saldat$month),
                      max(saldat$month)+60, 'year' ),
        col = rgb(0,0,0,0.1) )
points( value ~ month, data = saldat, pch = 21, lwd = 2, cex = 0.5,
        col = rgb(0,0,0,0.6), bg = rgb(1,1,1,1) )

# Plot frequencies
# select column
colnames( saldat )
var <- 'value'
x <- saldat[,var] |> na.omit()
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
obj |> plot( type = 'vectors', numvectors = 30 )
obj |> plot( type = 'paired', numvectors = 30 )
# W-correlation matrix
obj |> plot( type = 'wcor' )
wcor(obj,groups = 1:30) |> plot()


# SSA grouping
##
# Specify signal component groups
grp <- list( c(1,2,4,5), c(3,6,7), c(8,9), c(10,11,12,13), c(14,15,16,17) )
# grp <- list( c(1:7), c(8,9), c(10,11,12,13), c(14,15,16,17) )

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
ssa_sal_otb <- list( dat = data.frame( month = saldat$month,
                                        ts = saldat$value,
                                        signal = signal + mean(saldat$value),
                                        noise = noise + mean(saldat$value) ),
                      grp = grp,
                      sigstrength = sigstrength 
)
save( ssa_sal_otb, file = "../data/ssa_sal_otb.RData" )
