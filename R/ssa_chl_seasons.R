# OTB Pyrodinium signal processing with SSA
#
rm(list=ls(all=TRUE)) 

# Load libraries
if(!require(lubridate)) { install.packages('lubridate') }; library(lubridate)
if(!require(dplyr)) { install.packages('dplyr') }; library(dplyr)

# Load SSA outputs
load( "../data/ssa_chl_otb_nw.RData" )
load( "../data/ssa_chl_otb_ne.RData" )
load( "../data/ssa_chl_otb_cw.RData" )
load( "../data/ssa_chl_otb_ce.RData" )
load( "../data/ssa_chl_otb_sw.RData" )
load( "../data/ssa_chl_otb_se.RData" )

# Plot sub-segment signals by month
obj <- ls()[ c(4,3,2,1,6,5) ]
objlab <- c( "NW sub-segment", "NE sub-segment",
             "CW sub-segment", "CE sub-segment",
             "SW sub-segment", "SE sub-segment" )
colors <- hcl.colors( 24, "viridis", alpha = 0.7 )
png( "../figs/ssa_chl_seasons.png", height = 9, width = 9, units = 'in', res = 500 )
par( mfrow=c(3,2), mar=c(2,5,2,1) )
for( i in 1:length( obj ) ){
  # Subset sub-segment data
  this.i <- eval( parse(text=obj[i]) )$dat
  this.i$mo <- this.i$month |> month()
  # Initialize sub-segment plot
  plot( signal ~ mo, dat = this.i,
        type = 'l', ylim = c(0,30),
        las = 1, cex.axis = 1.3, cex.lab = 1.3,
        xlab = "", ylab = "Chlorophyll a (ug/L)", xaxt = 'n',
        lwd = 2, col = rgb(0,0,0,0) )
  axis( 1, at = 1:12, labels = month(1:12,label=TRUE), cex.axis = 1.3 )
  abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  mtext( objlab[i], side = 3, line = 0,
         font = 1, cex = 1, adj = 0 )
  # Loop through years to plot subseg signal by month
  yrs <- this.i$month |> year() |> unique() |> sort()
  for( j in 1:length(yrs) ){
    this.j <- this.i[ which( year(this.i$month)==yrs[j] ), ]
    lines( signal ~ mo, data = this.j,
           lwd = 3, col = colors[j] )
  }  # // end j loop
  # Draw legend in first plot
  if(i==1){
    polygon( x = c(1,2,2,1), y = c(30,30,20,20),
             col = rgb(1,1,1,1), border = rgb(0.3,0.3,0.3,1) )
    as.raster(colors,ncol=1) |> rasterImage(1,20,2,30)
    text( x = 2, y = 29.5, label = "2000", pos = 4 )
    text( x = 2, y = 20.5, label = "2023", pos = 4 )
  }
  
}  # // end i loop

dev.off()
