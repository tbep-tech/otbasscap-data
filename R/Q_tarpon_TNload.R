

rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/Q_tarpon.RData" )
load("../data-raw/totanndat.RData")

# Estimate TN loads for LTOC assuming TN conc is 1 mg/L (upr) and 0.5 mg/L (lwr)
tarpon$year <- year( tarpon$month )
tarpon_TN <- tarpon[,c('year','discharge')] |> group_by(year) |>
  dplyr::summarise(discharge_cfy=sum(discharge)) |> as.data.frame()
tarpon_TN$LTOC_TN_ton_upr <- tarpon_TN$discharge_cfy * 28.3168 * 1 * 1e-9  # 28.3 L/ft3, 1 mg TN/L, 1e-9 ton / mg
tarpon_TN$LTOC_TN_ton_lwr <- tarpon_TN$discharge_cfy * 28.3168 * 0.5 * 1e-9 # 28.3 L/ft3, 0.5 mg TN/L, 1e-9 ton / mg

# Generate plot of TN loads for LTOC
png( "../figs/Q_tarpon_TNload.png", width = 7, height = 5,
     units = 'in', res = 500 )
par( mar=c(2,4,1,1) )
plot( LTOC_TN_ton_upr ~ year, data = tarpon_TN,
      type ='l', lwd = 4, col = rgb(0,0.4,0.8,0.8),
      # main = "Lake Tarpon TN load per year (rough estimates)",
      ylab = "TN load (tons/year)", ylim = c(0,120), las = 1, xlab = "",
      cex.axis = 1.3, cex.lab = 1.3 )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
abline( h = seq(0,150,10), col = rgb(0,0,0,0.1) )
abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
abline( h = (40*86000*30*12) * 28.3168 * 1 * 1e-9,  # discharge 40 cfs & 1 mg/L TN
        col = rgb(0,0,0,0.6), lwd = 3 )
points( LTOC_TN_ton_upr ~ year, data = tarpon_TN, lwd = 3,
        pch = 21, bg = rgb(1,1,1,1), col = rgb(0,0.4,0.8,1) )
lines( LTOC_TN_ton_lwr ~ year, data = tarpon_TN,
       col = rgb(0,0.4,0.8,0.8), lty = 2, lwd = 2 )
points( LTOC_TN_ton_lwr ~ year, data = tarpon_TN, lwd = 2,
        pch = 21, bg = rgb(1,1,1,1), col = rgb(0,0.4,0.8,1) )
abline( h = (40*86000*30*12) * 28.3168 * 0.5 * 1e-9,  # discharge 40 cfs & 0.5 mg/L TN
        col = rgb(0,0,0,0.6), lty = 2, lwd = 2 )
legend( 'topright', bty = 'n', lwd = c(4,2), lty = c(1,2),
        col = c( rgb(0,0.4,0.8,0.8), rgb(0,0.4,0.8,0.8),
                 rgb(0,0,0,0.6), rgb(0,0,0,0.6) ),
        legend = c("Observed discharge & 1 mg/L TN",
                   "Observed discharge &  0.5 mg/L TN",
                   "40 cfs discharge & 1 mg/L TN",
                   "40 cfs discharge & 0.5 mg/L TN"
                   ) )
dev.off()

# Calculate TN load reductions associated with LTOC discharge reduction (40 cfs)
otb_TN <- totanndat[ which( totanndat$bay_segment=="Old Tampa Bay"),
                     c('year','tn_load') ] |> as.data.frame()
colnames(otb_TN)[2] <- 'OTB_TN_tons'
loads_TN <- left_join( tarpon_TN, otb_TN, 'year' )
loads_TN$OTB_TN_excess <- loads_TN$OTB_TN_tons - 486  # actual loads minus TMDL
loads_TN$LTOC_reduc_upr <- loads_TN$LTOC_TN_ton_upr - (40*86000*30*12)*28.3168*1*1e-9  # LTOC load reduction (40 cfs & 1 mg/L)
loads_TN$LTOC_reduc_lwr <- loads_TN$LTOC_TN_ton_lwr - (40*86000*30*12)*28.3168*0.5*1e-9  # LTOC load reduction (40 cfs & 0.5 mg/L)
loads_TN$LTOC_ratio_upr <- loads_TN$LTOC_reduc_upr / loads_TN$OTB_TN_excess # ratio of LTOC reduc to excess (40 cfs & 1 mg/L)
loads_TN$LTOC_ratio_lwr <- loads_TN$LTOC_reduc_lwr / loads_TN$OTB_TN_excess# ratio of LTOC reduc to excess (40 cfs & 0.5 mg/L)
loads_TN