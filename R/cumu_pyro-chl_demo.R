rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )
load( "../data/Pyro.Rdata")


png( "../figs/cumu_pyro-chl_demo.png", width = 7, height = 10, units = 'in', res = 500 )
par(mfrow=c(2,1))

# Assemble data for analysis
# chlorophyll data (EPC)
epcwq3.sub <- epcwq3[ which( epcwq3$param=="Chla" &
                               year(epcwq3$date) >= 2012 ), ]
epcwq3.sub$month <- floor_date( epcwq3.sub$date, unit = 'month' ) 
chldat <- epcwq3.sub |> group_by(month) |> summarise( chl = mean(value) ) |> as.data.frame()
# pyro data (FWC)
pyro.sub <- pyro[ which( pyro$yr>=2012 ), ]
pyro.sub$month <- floor_date( pyro.sub$date, unit = 'month' ) 
pyro.sub <- pyro.sub[ which(complete.cases(pyro.sub)), ]  # pyro==NA means zero cells/L
pyro.sub$logval <- log10( pyro.sub$pyro )
pyrodat <- pyro.sub |> group_by(month) |> summarise( pyro = max(logval) ) |> as.data.frame()
# join pyro and chl data by month
pcdat <- inner_join( pyrodat, chldat, by = 'month' )

# Plot histogram of max pyro assoc with mean chl-a at or below 9.3 ug/L
hist( pcdat$pyro[ which(pcdat$chl <= 9.3) ], breaks = 20, freq = TRUE,
      border = rgb(1,1,1,1), col = rgb(0,0,0,0.4),
      main = "A. Distribution of monthly max Pyro associated with\nmonthly mean chlorophyll-a at or below 9.3 ug/L",
      xlab = "cells/L", xaxt = 'n', xlim = c(1,7)
      )
axis( 1, at = 1:7, labels = c( expression(10^1), expression(10^2),
                               expression(10^3), expression(10^4),
                               expression(10^5), expression(10^6),
                               expression(10^7) ) )
abline( v = pcdat$pyro[ which(pcdat$chl <= 9.3) ] |> median(), 
        col = rgb(0,0.4,0.8,0.7), lwd = 2 )
abline( v = pcdat$pyro[ which(pcdat$chl <= 9.3) ] |> quantile(0.25),
        lty = 2, col = rgb(0,0.4,0.8,0.7), lwd = 2 )
abline( v = pcdat$pyro[ which(pcdat$chl <= 9.3) ] |> quantile(0.75),
        lty = 2, col = rgb(0,0.4,0.8,0.7), lwd = 2 )

# Define function to summarize pyro distribution
distn <- function( x, max_chl ){
  # subset pyro data associated with chl values at or below max_chl
  this <- x$pyro[ which( x$chl <= max_chl ) ]
  # assemble output statistics
  out <- data.frame( median = median(this),
                     min = min(this),
                     lwr_iqr = quantile(this, 0.25),
                     upr_iqr = quantile(this, 0.75),
                     max = max(this)
  )
  return( out )
}  # // end distn()


# Calculate pyro distn statistics
# initiate dataframe
pyro_chl <- data.frame( chl = seq(4,25,0.1),
                        median = NA,
                        min = NA,
                        lwr_iqr = NA,
                        upr_iqr = NA,
                        max = NA
)
# populate rows
for( i in 1:nrow(pyro_chl) ){
  pyro_chl[i,2:6] <- distn( pcdat, max_chl = pyro_chl[i,1] )
}  # // end i loop

# Plot pyro distn as a function of chl_max
plot( median ~ chl, data = pyro_chl, las = 1,
      type = 'l', lwd = 2, col = rgb(0,0.2,0.6,0.7),
      main = "B. OTB: Max Pyro distn ~ Mean Chl-a",
      ylim = c(1,7), yaxt = 'n',
      ylab = "Pyro (cells/L)", xlab = ""
)
mtext( "Chlorophyll a (ug/L)", side = 1, line = 2 )
axis( 2, at = 1:7, las = 1,
      labels = c( expression(10^1), expression(10^2),
                  expression(10^3), expression(10^4),
                  expression(10^5), expression(10^6),
                  expression(10^7) ) )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
abline( v = min(pyro_chl$chl):max(pyro_chl$chl), col = rgb(0,0,0,0.1) )
polygon( x = c( pyro_chl$chl, rev(pyro_chl$chl) ),
         y = c( pyro_chl$min, rev(pyro_chl$max) ),
         col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
polygon( x = c( pyro_chl$chl, rev(pyro_chl$chl) ),
         y = c( pyro_chl$lwr_iqr, rev(pyro_chl$upr_iqr) ),
         col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
segments( x0 = 9.3,
          y0 = 0, y1 = pyro_chl$upr_iqr[which(pyro_chl$chl==9.3)],
          lty = 2, col = 2 )
text( x = 9.3, y = pyro_chl$upr_iqr[which(pyro_chl$chl==9.3)],
      col = 2, labels = "9.3 ug/L", pos = 4, srt = 90 )
segments( x0 = 0, x1 = 9.3,
          y0 = pyro_chl$median[which(pyro_chl$chl==9.3)],
          lty = 2, col = 2 )
segments( x0 = 0, x1 = 9.3,
          y0 = pyro_chl$lwr_iqr[which(pyro_chl$chl==9.3)],
          lty = 2, col = 2 )
segments( x0 = 0, x1 = 9.3,
          y0 = pyro_chl$upr_iqr[which(pyro_chl$chl==9.3)],
          lty = 2, col = 2 )
legend( 'bottomright', bty = 'n',
        legend = c("Median of maxima",
                   "IQR of maxima",
                   "Range of maxima (min/max)"), text.font = 2,
        text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )


dev.off()