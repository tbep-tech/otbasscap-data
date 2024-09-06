# OTB Pyrodinium signal processing with SSA
#
rm(list=ls(all=TRUE)) 

# Load libraries
if(!require(lubridate)) { install.packages('lubridate') }; library(lubridate)
if(!require(dplyr)) { install.packages('dplyr') }; library(dplyr)

# Load SSA outputs
load( "../data/ssa_pyro_otb.RData" )
load( "../data/ssa_pyro_otb_nw.RData" )
load( "../data/ssa_pyro_otb_ne.RData" )
load( "../data/ssa_pyro_otb_cw.RData" )
load( "../data/ssa_pyro_otb_ce.RData" )
load( "../data/ssa_pyro_otb_sw.RData" )
load( "../data/ssa_pyro_otb_se.RData" )

# Plot signals and noise
obj <- ls()[ c(1,5,4,3,2,7,6) ]
objlab <- c( "OTB segment",
             "OTB, NW subsegment", "OTB, NE subsegment",
             "OTB, CW subsegment", "OTB, CE subsegment",
             "OTB, SW subsegment", "OTB, SE subsegment" )
png( "../figs/ssa_pyro_plots.png", height = 9, width = 7, units = 'in', res = 500 )
par( mfrow=c(7,1), mar=c(1,5,2,1), oma=c(2,0,0,0) )
for( i in 1:length( obj ) ){
  
  this <- eval( parse(text=obj[i]) )
  
  plot( ts ~ month, dat = this$dat,
        type = 'l', xlim = c( as.Date("2012-01-01"), as.Date("2023-12-31") ),
        ylim = c(0,7), las = 1,
        xlab = "", ylab = "", yaxt = 'n', xaxt = 'n',
        lwd = 2, col = rgb(0,0,0,0.6) )
  mtext( paste0( "(",letters[i], ") ", objlab[i],"  (signal strength: ",
                 round(this$sigstrength,2), "%)" ),
         side = 3, line = 0,
         font = 1, cex = 1.1, adj = 0 )
  date.seq <- seq.Date( min(this$dat$month), max(this$dat$month)+365, 'year' )
  if( i==length(obj) ){
    axis( 1, at = date.seq, labels = year(date.seq), cex.axis = 1.3  )
  } else {
    axis( 1, at = date.seq, labels = rep("",length(date.seq))  )
  }
  if( i==4 ){
    mtext( "cells/L",
           side = 2, line = 3, cex = 1.1 )
  }
  axis( 2, at = c(0,2,4,6), las = 1, cex.axis = 1.2,
        labels = c( 0, expression(10^2),
                    expression(10^4),
                    expression(10^6) ) )
  lines( signal ~ month, dat = this$dat,
         lwd = 4, col = rgb(1,0.2,0.1,0.8) )
  abline( v = date.seq, col = rgb(0,0,0,0.1))
  abline( h = axTicks(2), col = rgb(0,0,0,0.1))
  
  if( i==1 ){
    legend( "topright", horiz = TRUE, bty = 'n',
            xpd = TRUE, inset = c(0,-0.35),
            legend = c("data","signal"), lwd = c(2,4), cex = 1.1,
            col = c( rgb(0,0,0,0.6), rgb(1,0.2,0.1,0.8) ) )
  }
  
}  # // end i loop

dev.off()
