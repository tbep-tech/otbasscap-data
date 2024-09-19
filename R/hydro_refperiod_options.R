rm(list=ls(all=TRUE)) 

setwd("C:\\Users\\miles\\Documents\\GitHub\\otbasscap-data\\R")
load("../data-raw/totanndat.RData")

otbdat <- totanndat[ which( totanndat$bay_segment=="Old Tampa Bay"), ] |> as.data.frame()



# Compute the current hydrologically normalized TN load (Option 0)
hy_ref <- 449
otbdat$option0_tn_load_norm <- otbdat$tn_load * hy_ref / otbdat$hy_load


# Option 2. Updated TMDL and hydro reference period
ref_TN <- 462  # tons (annual average TN load 2005-2009)
ref_hydro <- 536.96  # million cubic meters (average annual hydro load 2005-2009)
otbdat$option2_tn_load_norm <- otbdat$tn_load * ref_hydro / otbdat$hy_load
png( "../figs/hydro_refperiod_option2.png", width = 8, height = 4, units = "in", res = 600 )
par( mfrow=c(1,1), mar=c(4,5,1,1) )
plot( tn_load ~ year, data = otbdat, type = 'l',
      main = "", xlab = "", ylab = "",
      ylim = c(0,1000),
      las = 1, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.3,
      col = rgb(0,0.1,0.4,0.2), lwd = 3 )
abline( v = axTicks(1), col = rgb(0,0,0,0.2) )  # grid
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )  # grid
points( tn_load ~ year, data = otbdat,  # actual load
        pch = 21, col = rgb(0,0.1,0.4,0.4),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )
mtext( "TN load (tons)", side = 2, line = 3.5, cex = 1.3 )  # ylab
abline( h = ref_TN, col = rgb(0.4,0.7,0.4,0.5), lty = 2, lwd = 3 )  # updated TMDL
abline( h = 486, col = rgb(0.7,0.4,0.4,0.5), lty = 3, lwd = 2 )  # current TMDL
lines( option0_tn_load_norm  ~ year, data = otbdat,  # current normalized load
       col = rgb(0.7,0.1,0.2,0.8), lwd = 3 )
points( option0_tn_load_norm  ~ year, data = otbdat, pch = 21, col = rgb(0.7,0.1,0.2,0.8),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

lines( option2_tn_load_norm  ~ year, data = otbdat,  # Option 2 normalized load
       col = rgb(0.1,0.7,0.2,0.8), lwd = 3 )
points( option2_tn_load_norm  ~ year, data = otbdat, pch = 21, col = rgb(0.1,0.7,0.2,0.8),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

legend( 'bottomleft', horiz = FALSE, bty='n',
        lty = c(1,1,3,1,2), lwd = c(3,3,2,3,3), cex = 1,
        col = c( rgb(0,0.1,0.4,0.2), rgb(0.7,0.1,0.2,0.8), rgb(0.7,0.4,0.4,0.5),
                 rgb(0.1,0.7,0.2,0.8), rgb(0.4,0.7,0.4,0.5) ),
        legend = c("Actual load",
                   "Normalized load (1992-1994 ref period)",
                   expression("Current TMDL (486 tons/yr)"),
                   "Normalized load (2005-2009 ref period)",
                   expression("Updated TMDL (462 tons/yr)") )
                    
)
dev.off()


# Option 3. Current TMDL and updated hydro reference period
otbdat$option3_tn_load_norm <- otbdat$tn_load * ref_hydro / otbdat$hy_load
png( "../figs/hydro_refperiod_option3.png", width = 8, height = 4, units = "in", res = 600 )
par( mfrow=c(1,1), mar=c(4,5,1,1) )
plot( tn_load ~ year, data = otbdat, type = 'l',
      main = "", xlab = "", ylab = "",
      ylim = c(0,1000),
      las = 1, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.3,
      col = rgb(0,0.1,0.4,0.2), lwd = 3 )
abline( v = axTicks(1), col = rgb(0,0,0,0.2) )  # grid
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )  # grid
points( tn_load ~ year, data = otbdat,  # actual load
        pch = 21, col = rgb(0,0.1,0.4,0.4),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )
mtext( "TN load (tons)", side = 2, line = 3.5, cex = 1.3 )  # ylab
abline( h = 486, col = rgb(0,0,0,0.5), lty = 2, lwd = 3 )  # TMDL

lines( option0_tn_load_norm  ~ year, data = otbdat,  # current normalized load
       col = rgb(0.7,0.1,0.2,0.8), lwd = 3 )
points( option0_tn_load_norm  ~ year, data = otbdat, pch = 21, col = rgb(0.7,0.1,0.2,0.8),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

lines( option3_tn_load_norm  ~ year, data = otbdat,  # Option 3 normalized load
       col = rgb(0.1,0.7,0.2,0.8), lwd = 3 )
points( option3_tn_load_norm  ~ year, data = otbdat, pch = 21, col = rgb(0.1,0.7,0.2,0.8),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

legend( 'bottomleft', horiz = FALSE, bty='n',
        lty = c(1,1,1,2), lwd = 3, cex = 1,
        col = c( rgb(0,0.1,0.4,0.2), rgb(0.7,0.1,0.2,0.8),
                 rgb(0.1,0.7,0.2,0.8), rgb(0,0,0,0.5) ),
        legend = c("Actual load","Normalized load (1992-1994 ref period)",
                   "Normalized load (1992-1994 TN ref period; 2005-2009 hydro ref period)",
                   expression("TMDL (486 tons/yr)") )
)
dev.off()

# Option 4
TN_safety <- 522  # 522 tons/year is associated with 50% lower bound attainment probability, P(chl<8.5)
png( "../figs/hydro_refperiod_option4.png", width = 8, height = 4, units = "in", res = 600 )
par( mfrow=c(1,1), mar=c(4,5,1,1) )
plot( tn_load ~ year, data = otbdat, type = 'l',
      main = "", xlab = "", ylab = "",
      ylim = c(0,1000),
      las = 1, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.3,
      col = rgb(0,0.1,0.4,0.8), lwd = 3 )
abline( v = axTicks(1), col = rgb(0,0,0,0.2) )  # grid
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )  # grid
mtext( "TN load (tons)", side = 2, line = 3.5, cex = 1.3 )  # ylab
abline( h = TN_safety, col = rgb(0,0.7,1,0.5), lty = 2, lwd = 3 )  # TMDL+MOS
abline( h = 486, col = rgb(0,0,0,0.5), lty = 2, lwd = 3 )  # TMDL
lines( option0_tn_load_norm  ~ year, data = otbdat,  # current normalized load
       col = rgb(0.8,0.1,0.2,0.2), lwd = 3 )
points( option0_tn_load_norm  ~ year, data = otbdat,
        pch = 21, col = rgb(0.8,0.1,0.2,0.4),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )
lines( tn_load ~ year, data = otbdat,  # actual load
       col = rgb(0,0.1,0.4,0.8), lwd = 3 )
points( tn_load ~ year, data = otbdat,  # actual load
        pch = 21, col = rgb(0,0.1,0.4,0.8),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

legend( 'bottomleft', horiz = FALSE, bty='n',
        lty = c(1,1,2,2), lwd = 3, cex = 1,
        col = c( rgb(0,0.1,0.4,0.8), rgb(0.8,0.1,0.2,0.4),
                 rgb(0,0.7,1,0.5), rgb(0,0,0,0.5) ),
        legend = c("Actual load",
                   "Normalized load (1992-1994 ref period)",
                   expression("TMDL + safety factor (522 tons/yr)"),
                   expression("TMDL (486 tons/yr)") )
)
dev.off()