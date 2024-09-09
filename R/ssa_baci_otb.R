# OTB Bacillaria signal processing with SSA
#
rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load & process phytoplankton data
load( "../data-clean/epcphyto.Rdata" )
# Subset Bacillaria data
phyto.sub <- phyto[ which( year(phyto$date) >= 2012 &
                             year(phyto$date) <= 2023 ), ]
phyto.sub <- phyto.sub[ which( phyto.sub$family=="Bacillariaceae" ), ]
# Aggregate cell counts to monthly scale
phyto.sub$month <- floor_date( phyto.sub$date, unit = 'month' )
phytodat <- phyto.sub |> group_by(month) |>
  dplyr::summarise( cells = max(cellcount) ) |> as.data.frame()
date.seq <- data.frame( month = seq.Date(as.Date("2012-01-01"),
                                         as.Date("2023-01-01"),'month') )
baci <- left_join( date.seq, phytodat, 'month' )
# interpolate missing values
missing.idx <- baci$cells |> is.na() |> which()
for( i in 1:length(missing.idx) ){  # this loop is specific to these data and..
  if( !is.na(baci$cells[missing.idx[i]+1]) ){  # ..does not handle all possible cases
    baci$cells[missing.idx[i]] <- mean( c(baci$cells[missing.idx[i]-1],
                                          baci$cells[missing.idx[i]+1]) )
  } else {
    baci$cells[missing.idx[i]] <- baci$cells[missing.idx[i]-1] +
      diff(c( baci$cells[missing.idx[i]-1],
              baci$cells[missing.idx[i]+2] ))/3
    baci$cells[missing.idx[i]+1] <- baci$cells[missing.idx[i]+2] -
      diff(c( baci$cells[missing.idx[i]-1],
              baci$cells[missing.idx[i]+2] ))/3
  }
}  # // end interpolation loop
baci$cells |> is.na() |> any()  # confirm FALSE
# Log-transform cell counts
baci$logcells <- log10( baci$cells + 1 )


# Plot time series
plot( logcells ~ month, data = baci, las = 1,
      type = 'l', lwd = 2, col = rgb(0,0,0,0.6),
      main = "Bacillaria ts", xlab = '', ylab = 'Bacillaria (log10 cells/L)' )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
abline( v = seq.Date( min(baci$month),
                      max(baci$month)+60, 'year' ),
        col = rgb(0,0,0,0.1) )
points( logcells ~ month, data = baci, pch = 21, lwd = 2, cex = 0.5,
        col = rgb(0,0,0,0.6), bg = rgb(1,1,1,1) )
points( logcells ~ month, data = baci[missing.idx,], col = 2, pch = 16 )
legend( 'topleft', bty='n', legend='interpolated value', pch=16, col=2 )


# Plot frequencies
# select column
colnames( baci )
var <- 'logcells'
x <- baci[,var] |> na.omit()
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
win <- 60
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
wcor(obj,groups = 1:20) |> plot()


# SSA grouping
##
# Specify signal component groups
grp <- list( c(1,2), c(3,4,6,7) )

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
ssa_baci_otb <- list( dat = data.frame( month = baci$month,
                                        ts = baci$logcells,
                                        signal = signal + mean(baci$logcells),
                                        noise = noise + mean(baci$logcells) ),
                                        grp = grp,
                                        sigstrength = sigstrength 
)
save( ssa_baci_otb, file = "../data/ssa_baci_otb.RData" )  
  
