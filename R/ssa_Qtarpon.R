# OTB Chlorophyll-a signal processing with SSA
#
rm(list=ls(all=TRUE)) 

# Load libraries
if(!require(Rssa)) { install.packages('Rssa') }; library(Rssa)
if(!require(lubridate)) { install.packages('lubridate') }; library(lubridate)
if(!require(dplyr)) { install.packages('dplyr') }; library(dplyr)

# Load & process data
load( "../data-clean/Q_tarpon.Rdata")

# Plot time series
par(mfrow=c(2,1), mar=c(4,5,4,1))
plot( discharge ~ month, data = tarpon, las = 1,
      type = 'l', lwd = 2, col = rgb(0,0,0,0.6), ylim = c(0,max(tarpon$discharge)),
      main = "Lake Tarpon canal discharge ts", xlab = '', ylab = 'ft3/mo' )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
abline( v = seq.Date( min(tarpon$month),
                      max(tarpon$month)+60, 'year' ),
        col = rgb(0,0,0,0.1) )
points( discharge ~ month, data = tarpon, pch = 21, lwd = 2, cex = 0.5,
        col = rgb(0,0,0,0.6), bg = rgb(1,1,1,1) )

# Plot frequencies
# select column
colnames( tarpon )
var <- 'discharge'
x <- tarpon[,var] |> na.omit()
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
win <- 132
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
grp <- list( c(1,2,6,10), c(3,4) )

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
ssa_Qtarpon_otb <- list( dat = data.frame( month = tarpon$month,
                                           ts = tarpon$discharge,
                                           signal = signal + mean(tarpon$discharge),
                                           noise = noise + mean(tarpon$discharge) ),
                         grp = grp,
                         sigstrength = sigstrength 
)
save( ssa_Qtarpon_otb, file = "../data/ssa_Qtarpon_otb.RData" )


# Generate figure
png( "../figs/ssa_Qtarpon_otb.png",
     width = 7, height = 10, res = 500, units = "in" )
par( mfrow=c(4,1), mar=c(0,5,2,1), oma=c(2,0,0,0) )
# data and signal
plot( ts ~ month, dat = ssa_Qtarpon_otb$dat,
      type = 'l', las = 1,
      xlab = "", xaxt = 'n',
      ylab = "", cex.lab = 1.3,
      ylim = c(0,max(ssa_Qtarpon_otb$dat$ts)), cex.axis = 1.3,
      lwd = 2, col = rgb(0,0,0,0.8) )
mtext( paste0( "(a) Discharge at Lake Tarpon canal  (signal strength: ",
               round(ssa_Qtarpon_otb$sigstrength,2), "%)" ),
       side = 3, line = 0,
       font = 1, cex = 1.1, adj = 0 )
date.seq <- seq.Date( min(ssa_Qtarpon_otb$dat$month),
                      max(ssa_Qtarpon_otb$dat$month)+365, 'year' )
axis( 1, at = date.seq,
      labels = rep("",length(date.seq)), cex.axis = 1.3  )
lines( signal ~ month, dat = ssa_Qtarpon_otb$dat,
       lwd = 4, col = rgb(1,0.2,0.1,0.8) )
abline( v = date.seq, col = rgb(0,0,0,0.1))
abline( h = axTicks(2), col = rgb(0,0,0,0.1))
legend( "topright", horiz = TRUE, bty = 'n',
        legend = c("data","signal"), lwd = c(2,4), cex = 1.3,
        col = c( rgb(0,0,0,0.8), rgb(1,0.2,0.1,0.8) ) )
# signal components
ylims <- range(unlist(recon))
for( i in 1:length(recon) ){
  if( i < length(grp) ){
    # Plot signal components
    plot( recon[[i]] ~ ssa_Qtarpon_otb$dat$month, type = 'l',
          col = rgb(1,0.2,0.1,0.8), lwd = 2,
          ylim = ylims, las = 1, cex.axis = 1.3,
          main = "", xlab = '', ylab = '', xaxt = 'n'
    )
    if(i==1){
      mtext( "(b) Discharge time series components",
             side = 3, line = 0,
             font = 1, cex = 1.1, adj = 0 )
    }
    if(i==2){
      mtext( "Discharge (ft3/mo)",
             side = 2, line = 3,
             font = 1, cex = 1.1 )
    }
    abline( v = date.seq, col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    rspec <- spectrum(recon[[i]], method = 'pgram', plot = FALSE )
    rspec <- data.frame( power = rspec$spec, period = 1/rspec$freq )
    rspec <- rspec[ order( rspec$power, decreasing = TRUE ), ]
    peak <- rspec$period[1] |> round(1)
    legend( 'topright', bty = 'n', cex = 1.3,
            legend = paste0("Signal component ",i,
                            "  (",round(varexp[[i]],2),"%)",
                            "\nPeak wavelength: ",peak," months"))
    axis( 1, at = date.seq,
          labels = rep("",length(date.seq)), cex.axis = 1.3  )
  } else {
    # Plot noise
    plot( recon[[i]] ~ ssa_Qtarpon_otb$dat$month, type = 'l',
          col = rgb(0.3,0.3,0.3,0.8), lwd = 2, cex.axis = 1.3,
          ylim = ylims, las = 1, xaxt = 'n',
          main = "",
          xlab = '', ylab = ''
    )
    abline( v = date.seq, col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    legend( 'topright', bty = 'n', cex = 1.3,
            legend = paste0("Noise","  (",round(varexp[[i]],2),"%)") )
    axis( 1, at = date.seq,
          labels = year(date.seq), cex.axis = 1.3  )
  }
}  # // end i 
dev.off()
