

rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/Q_tarpon.RData" )

# Estimate TN loads assuming TN conc is 1 mg/L and 0.5 mg/L
tarpon$year <- year( tarpon$month )
loads <- tarpon[,c('year','discharge')] |> group_by(year) |>
  dplyr::summarise(discharge_cfy=sum(discharge)) |> as.data.frame()
loads$TN_ton_upr <- loads$discharge_cfy * 28.3168 * 1 * 1e-9  # 28.3 L/ft3, 1 mg TN/L, 1e-9 ton / mg
loads$TN_ton_lwr <- loads$discharge_cfy * 28.3168 * 0.5 * 1e-9 # 28.3 L/ft3, 0.5 mg TN/L, 1e-9 ton / mg

# Plot
png( "../figs/Q_tarpon_TNload.png", width = 7, height = 5,
     units = 'in', res = 500 )
par( mar=c(2,4,1,1) )
plot( TN_ton_upr ~ year, data = loads,
      type ='l', lwd = 4, col = rgb(0,0.4,0.8,0.8),
      # main = "Lake Tarpon TN load per year (rough estimates)",
      ylab = "TN load (tons/year)", ylim = c(0,120), las = 1, xlab = "",
      cex.axis = 1.3, cex.lab = 1.3 )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
abline( h = seq(0,150,10), col = rgb(0,0,0,0.1) )
abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
abline( h = (40*86000*30*12) * 28.3168 * 1 * 1e-9,  # discharge 40 cfs & 1 mg/L TN
        col = rgb(0,0,0,0.6), lwd = 3 )
points( TN_ton_upr ~ year, data = loads, lwd = 3,
        pch = 21, bg = rgb(1,1,1,1), col = rgb(0,0.4,0.8,1) )
lines( TN_ton_lwr ~ year, data = loads,
       col = rgb(0,0.4,0.8,0.8), lty = 2, lwd = 2 )
points( TN_ton_lwr ~ year, data = loads, lwd = 2,
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
