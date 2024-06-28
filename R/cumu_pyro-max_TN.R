rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data/loads.RData" )
load( "../data/Pyro.Rdata")

# Specify subsegment
# subsegs <- c("NW","NE","CW","CE","SW","SE")

# Loop over sub-segments to generate plots
png( "../figs/cumu_pyro-max_TN.png", width = 7, height = 5, units = 'in', res = 500 )
# par( mfrow=c(3,4), mar=c(4,4,3,2) )
# for( subseg in subsegs ){
  
  # Assemble data for analysis
  # TN load data (TBEP)
  loaddat <- loads[ which(loads$param=="TN load"), ]
  loaddat <- loads[ which( year(loaddat$date) >= 2012 ), ]
  loaddat <- select( loaddat, date, value )
  colnames(loaddat) <- c("month","TN")
  
  # # sum loads over 3- month window
  # loaddat$TN_d1 <- c( loaddat$TN[2:nrow(loaddat)] , NA )
  # loaddat$TN_d2 <- c( loaddat$TN[3:nrow(loaddat)] , NA, NA )
  # loaddat$TN <- loaddat$TN + loaddat$TN_d1 + loaddat$TN_d2
  # loaddat <- loaddat[ which(complete.cases(loaddat)), ]
  # loaddat <- select( loaddat, month, TN )
  
  # pyro data (FWC)
  pyro.sub <- pyro[ which( pyro$yr>=2012
                           # & pyro$subsegment==subseg
                           ), ]
  pyro.sub$month <- floor_date( pyro.sub$date, unit = 'month' ) 
  pyro.sub <- pyro.sub[ which(complete.cases(pyro.sub)), ]  # pyro==NA means zero cells/L
  pyro.sub$logval <- log10( pyro.sub$pyro )
  pyrodat <- pyro.sub |> group_by(month) |> summarise( pyro = max(logval) ) |> as.data.frame()
  # join pyro and chl data by month
  pcdat <- inner_join( pyrodat, loaddat, by = 'month' )
  
  # Define function to summarize pyro distribution
  distn <- function( x, max_TN ){
    # subset pyro data associated with chl values at or below max_TN
    this <- x$pyro[ which( x$TN <= max_TN ) ]
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
  pyro_TN <- data.frame( TN = seq(10,300,10),
                         median = NA,
                         min = NA,
                         lwr_iqr = NA,
                         upr_iqr = NA,
                         max = NA
  )
  # populate rows
  for( i in 1:nrow(pyro_TN) ){
    pyro_TN[i,2:6] <- distn( pcdat, max_TN = pyro_TN[i,1] )
  }  # // end i loop
  
  # Plot pyro distn as a function of chl_max
  plot( median ~ TN, data = pyro_TN, las = 1,
        type = 'l', lwd = 2, col = rgb(0,0.2,0.6,0.7),
        main = paste0( 
                       # subseg,
                       " OTB\nMax Pyro distn ~ TN load" ),
        ylim = c(1,7), yaxt = 'n',
        ylab = "Pyro (cells/L)", xlab = ""
  )
  mtext( "TN load (tons/month)", side = 1, line = 2, cex = 0.7 )
  axis( 2, at = 1:7, las = 1,
        labels = c( expression(10^1), expression(10^2),
                    expression(10^3), expression(10^4),
                    expression(10^5), expression(10^6),
                    expression(10^7) ) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  abline( h = log10( c(2e0,3e0,4e0,5e0,6e0,7e0,8e0,9e0,
                       2e1,3e1,4e1,5e1,6e1,7e1,8e1,9e1,
                       2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,
                       2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,
                       2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,
                       2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,
                       2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6) ), col = rgb(0,0,0,0.1) )
  abline( v = pyro_TN$TN, col = rgb(0,0,0,0.1) )
  polygon( x = c( pyro_TN$TN, rev(pyro_TN$TN) ),
           y = c( pyro_TN$min, rev(pyro_TN$max) ),
           col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
  polygon( x = c( pyro_TN$TN, rev(pyro_TN$TN) ),
           y = c( pyro_TN$lwr_iqr, rev(pyro_TN$upr_iqr) ),
           col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
  segments( x0 = 40,
            y0 = 0, y1 = pyro_TN$upr_iqr[which(pyro_TN$TN==40)],
            lty = 2, col = 2 )
  text( x = 40, y = pyro_TN$upr_iqr[which(pyro_TN$TN==40)],
        col = 2, labels = "40 tons", pos = 4, srt = 90 )
  segments( x0 = 0, x1 = 40,
            y0 = pyro_TN$median[which(pyro_TN$TN==40)],
            lty = 2, col = 2 )
  segments( x0 = 0, x1 = 40,
            y0 = pyro_TN$lwr_iqr[which(pyro_TN$TN==40)],
            lty = 2, col = 2 )
  segments( x0 = 0, x1 = 40,
            y0 = pyro_TN$upr_iqr[which(pyro_TN$TN==40)],
            lty = 2, col = 2 )
  legend( 'bottomright', bty = 'n',
          legend = c("Median of maxima",
                     "IQR of maxima",
                     "Range of maxima (min/max)"), text.font = 2,
          text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )
  
  
# }  # end subseg loop

dev.off()
