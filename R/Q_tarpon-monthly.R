rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/Q_tarpon.RData" )

# Plot monthly discharge data
png( "../figs/Q_tarpon-monthly.png", width = 7, height = 3,
     units = 'in', res = 500 )
par( mar=c(2,4,1,1) )
  plot( discharge/1e6 ~ month, data = tarpon,
        type = 'l', las = 1,
        col = rgb(0,0.4,0.8,0.8), lwd = 2,
        main = "Monthly Discharge from Lake Tarpon",
        ylab = "Discharge (million ft3/mo)", xlab = '' )
  abline( v = seq.Date( floor_date(min(tarpon$month),'year'),
                        ceiling_date(max(tarpon$month),'year'), 'year' ),
          col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  points( discharge/1e6 ~ month, data = tarpon, cex = 0.4,
          pch = 21, bg = rgb(1,1,1,1), col = rgb(0,0.4,0.8,1) )
dev.off()
