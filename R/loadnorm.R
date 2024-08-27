rm(list=ls(all=TRUE)) 

setwd("C:\\Users\\miles\\Documents\\GitHub\\otbasscap-data\\R")
load("../data-raw/totanndat.RData")

otbdat <- totanndat[ which( totanndat$bay_segment=="Old Tampa Bay"), ] |> as.data.frame()
head( otbdat )



# Compute hydrologically normalized TN load
hy_ref <- 449
otbdat$tn_load_norm <- otbdat$tn_load * hy_ref / otbdat$hy_load

# Plot TN loads
png( "../figs/loadnorm.png", width = 8, height = 8, units = "in", res = 600 )
par( mfrow=c(2,1), mar=c(4,5,3,5) )
plot( tn_load ~ year, data = otbdat, type = 'l',
      main = "(a) TN loads to OTB", xlab = "", ylab = "",
      ylim = c(0,1000),
      las = 1, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.3,
      col = rgb(0,0.1,0.4,0.8), lwd = 3 )

abline( v = axTicks(1), col = rgb(0,0,0,0.2) )  # grid
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )  # grid

points( tn_load ~ year, data = otbdat,  # actual load
        pch = 21, col = rgb(0,0.1,0.4,1),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

mtext( "TN load (tons)", side = 2, line = 3.5, cex = 1.3 )  # ylab
mtext( expression("TN load rate (tons /"~10^6~"m3)"),
       side = 4, line = 3.5, cex = 1.3 )  # y2 lab
axis( 4, at = seq(0,1000,200), labels = round(seq(0,2.2,length.out=6),1),     # y2
      cex.lab = 1.3, cex.axis = 1.3, las = 1)

abline( h = 486, col = rgb(0,0,0,0.5), lty = 2, lwd = 3 )  # TMDL

lines( tn_load_norm ~ year, data = otbdat,  # normalized load
       col = rgb(0.7,0.1,0.2,0.8), lwd = 3 )
points( tn_load_norm ~ year, data = otbdat, pch = 21, col = rgb(0.7,0.1,0.2,0.8),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

legend( 'bottomleft', horiz = FALSE, bty='n',
        lty = c(1,1,2), lwd = 3, cex = 1.2,
        col = c( rgb(0,0.1,0.4,0.8), rgb(0.7,0.1,0.2,0.8), rgb(0,0,0,0.5) ),
        legend = c("Actual load","Normalized load",
                   expression("TMDL (486 tons/yr) = Delivery ratio (1.08 tons /"*10^6~"m3)") )
         )


# Plot hydro loads
plot( hy_load ~ year, data = otbdat, type = 'l',
      main = "(b) Hydrologic loads to OTB", xlab = "", ylab = "",
      las = 1, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.3,
      col = rgb(0,0.4,0.7,0.8), lwd = 3 )

abline( v = axTicks(1), col = rgb(0,0,0,0.2) )  # grid
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )  # grid

abline( h = hy_ref, col = rgb(0,0,0,0.5), lty = 2, lwd = 3 )  # reference hydrologic load
text( x = 2012, y = hy_ref, pos = 1, col = rgb(0,0,0,0.7),
      labels = "Reference hydrologic load (449 million m3)" )

points( hy_load ~ year, data = otbdat, pch = 21, col = rgb(0,0.4,0.7,1),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

mtext( "million m3", side = 2, line = 3.5, cex = 1.3 )  # ylab



dev.off()