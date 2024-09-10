# Cond dist plot for chl~TN during PYRO BLOOM months only

rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )
load( "../data/loads.Rdata")
load( "../data/Pyro.Rdata")

# Get non-zero Pyro months
pyro.sub <- pyro[ which( pyro$yr>=2012 ), ]
pyro.sub$month <- floor_date( pyro.sub$date, unit = 'month' ) 
pyro.months <- pyro.sub$month[ which(complete.cases(pyro.sub)) ] |>
  unique() |> as.data.frame()  # pyro==NA means zero cells/L
colnames(pyro.months) <- "month"

# Specify subsegment
subsegs <- c("NW","NE","CW","CE","SW","SE")

# Loop over sub-segments to generate plots
png( "../figs/cumu_chl-mean_TN_otb-subsegments2.png",
     width = 7, height = 9, units = 'in', res = 500 )
par( mfrow=c(3,2), mar=c(4,4,2,2) )
for( subseg in subsegs ){

# Assemble data for analysis
# TN load data (TBEP)
loaddat <- loads[ which(loads$param=="TN load"), ]
loaddat <- select( loaddat, date, value )
colnames(loaddat) <- c("month","TN")

# # sum loads over 3- month window
# loaddat$TN_d1 <- c( loaddat$TN[2:nrow(loaddat)] , NA )
# loaddat$TN_d2 <- c( loaddat$TN[3:nrow(loaddat)] , NA, NA )
# loaddat$TN <- loaddat$TN + loaddat$TN_d1 + loaddat$TN_d2
# loaddat <- loaddat[ which(complete.cases(loaddat)), ]
# loaddat <- select( loaddat, month, TN )

# chlorophyll data (EPC)
epcwq3.sub <- epcwq3[ which( epcwq3$param=="Chla"
                             & epcwq3$subseg==subseg
                           ), ]
epcwq3.sub$month <- floor_date( epcwq3.sub$date, unit = 'month' ) 
chldat <- epcwq3.sub |> group_by(month) |> dplyr::summarise( chl = mean(value) ) |> as.data.frame()
# chldat$chl <- chldat$chl |> log10()
wqdat <- inner_join( pyro.months, chldat, by = 'month' )
wqdat <- inner_join( wqdat, loaddat, by = 'month' )

# Define function to summarize chl distribution
distn <- function( x, max_TN ){
  # subset chl data associated with chl values at or below max_TN
  this <- x$chl[ which( x$TN <= max_TN ) ]
  out <- data.frame( mean = mean(this),
                     min = min(this),
                     lwr = mean(this) - sd(this),
                     upr = mean(this) + sd(this),
                     max = max(this)
  )
  return( out )
}  # // end distn()


# Calculate chl distn statistics
# initiate dataframe
chl_TN <- data.frame( TN = seq(10,200,5),
                      mean = NA,
                      min = NA,
                      lwr = NA,
                      upr = NA,
                      max = NA
)
# populate rows
for( i in 1:nrow(chl_TN) ){
  chl_TN[i,2:6] <- distn( wqdat, max_TN = chl_TN[i,1] )
}  # // end i loop

# Plot chl distn as a function of chl_max
plot( mean ~ TN, data = chl_TN, las = 1,
      type = 'l', lwd = 2, col = rgb(0,0.2,0.6,0.7),
      main = paste0( 
        subseg," sub-segment"
        # " OTB\nMonthly mean Chl-a distn ~ TN load"
        ),
      ylim = c(0,50), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5,
      ylab = "", xlab = "", xlim = c(10,200)
)
if(subseg %in% c("SW","SE")){
  mtext( "TN load (tons/month)", side = 1, line = 2.5, cex = 1 )
}
if(subseg %in% c("CW")){
  mtext( "Chlorophyll a (ug/L)", side = 2, line = 2.5, cex = 1 )
}
# axis( 2, at = 0:3, las = 1,
      # labels = c( 10^0, 10^1, 10^2, 10^3 ) )
abline( h = seq(0,60,5), col = rgb(0,0,0,0.1) )
abline( h = seq(0,60,10), col = rgb(0,0,0,0.1) )
# abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
# abline( h = log10( c(1e0,2e0,3e0,4e0,5e0,6e0,7e0,8e0,9e0,
#                      1e1,2e1,3e1,4e1,5e1,6e1,7e1,8e1,9e1,
#                      1e2,2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2 ) ), col = rgb(0,0,0,0.1) )
abline( v = chl_TN$TN, col = rgb(0,0,0,0.1) )
polygon( x = c( chl_TN$TN, rev(chl_TN$TN) ),
         y = c( chl_TN$min, rev(chl_TN$max) ),
         col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
polygon( x = c( chl_TN$TN, rev(chl_TN$TN) ),
         y = c( chl_TN$lwr, rev(chl_TN$upr) ),
         col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )

segments( x0 = 40,
          y0 = -10, y1 = chl_TN$upr[which(chl_TN$TN==40)],
          lty = 2, col = 1 )
segments( x0 = 0, x1 = 40,
          y0 = chl_TN$mean[which(chl_TN$TN==40)],
          lty = 2, col = 1 )
segments( x0 = 0, x1 = 40,
          y0 = chl_TN$lwr[which(chl_TN$TN==40)],
          lty = 2, col = 1 )
segments( x0 = 0, x1 = 40,
          y0 = chl_TN$upr[which(chl_TN$TN==40)],
          lty = 2, col = 1 )

segments( x0 = 65,
          y0 = -10, y1 = chl_TN$upr[which(chl_TN$TN==65)],
          lty = 2, col = 2 )
segments( x0 = 0, x1 = 65,
          y0 = chl_TN$mean[which(chl_TN$TN==65)],
          lty = 2, col = 2 )
segments( x0 = 0, x1 = 65,
          y0 = chl_TN$lwr[which(chl_TN$TN==65)],
          lty = 2, col = 2 )
segments( x0 = 0, x1 = 65,
          y0 = chl_TN$upr[which(chl_TN$TN==65)],
          lty = 2, col = 2 )

text( x = 40, y = chl_TN$upr[which(chl_TN$TN==40)],
      col = 1, labels = "40 tons", pos = 4, srt = 90 )
text( x = 65, y = chl_TN$upr[which(chl_TN$TN==65)],
      col = 2, labels = "65 tons", pos = 4, srt = 90 )

legend( 'topleft', bty = 'n', cex = 1.3,
        legend = c("Mean of means",
                   "Mean\U00B1SD of means",
                   "Range of means (min/max)"), text.font = 2,
        text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )


}  # end subseg loop

dev.off()
