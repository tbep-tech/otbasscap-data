rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/Q_tarpon.RData" )
load( "../data/Pyro.Rdata")

# Subset routine Pyro samples
pyro <- pyro[ which(pyro$routine==TRUE), ]

# Overwrite 'discharge' column with log-transformed values
tarpon$discharge <- tarpon$logdischarge

# Specify subsegment
subsegs <- c("NW","NE","CW","CE","SW","SE")

# Loop over sub-segments to generate plots
png( "../figs/cumu_pyro-max_Qtarpon.png", width = 7, height = 9, units = 'in', res = 500 )
par( mfrow=c(3,2), mar=c(4,4,3,2) )
for( subseg in subsegs ){
  
  # Assemble data for analysis
  pyro.sub <- pyro[ which( pyro$subsegment==subseg & pyro$yr>=2012 ), ]
  pyro.sub$month <- floor_date( pyro.sub$date, unit = 'month' ) 
  pyro.sub <- pyro.sub[ which(complete.cases(pyro.sub)), ]  # pyro==NA means zero cells/L
  pyro.sub$logval <- log10( pyro.sub$pyro )
  pyrodat <- pyro.sub |> group_by(month) |> dplyr::summarise( pyro = max(logval) ) |> as.data.frame()
  # join pyro and discharge data by month
  pcdat <- inner_join( pyrodat, tarpon[,c('month','discharge')], by = 'month' )
  
  # Define function to summarize pyro distribution
  distn <- function( x, max_q ){
    # subset pyro data associated with discharge values at or below max_q
    this <- x$pyro[ which( x$discharge <= max_q ) ]
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
  pyro_q <- data.frame( discharge = log10( c( seq(1e1,9e1,1e1),
                                              seq(1e2,9e2,1e2),
                                              seq(1e3,9e3,1e3),
                                              seq(1e4,9e4,1e4),
                                              seq(1e5,9e5,1e5),
                                              seq(1e6,9e6,1e6),
                                              seq(1e7,9e7,1e7),
                                              seq(1e8,9e8,1e8),
                                              seq(1e9,9e9,1e9) ) ),
                        median = NA,
                        min = NA,
                        lwr_iqr = NA,
                        upr_iqr = NA,
                        max = NA
  )
  # populate rows
  for( i in 1:nrow(pyro_q) ){
    pyro_q[i,2:6] <- distn( pcdat, max_q = pyro_q[i,1] )
  }  # // end i loop
  
  # Plot pyro distn as a function of q_max
  plot( median ~ discharge, data = pyro_q, las = 1,
        type = 'l', lwd = 2, col = rgb(0,0.2,0.6,0.7),
        ylim = c(1,7), yaxt = 'n', xaxt = 'n',
        ylab = expression(italic(P.~bahamense)*" (cells/L)"), xlab = ""
  )
  mtext( paste0(subseg, " sub-segment"),
         side = 3, adj = 0, line = 1 )
  mtext( "Discharge (ft3/month)", side = 1, line = 2, cex = 0.7 )
  axis( 1, at = 1:9, las = 1,
        labels = c( expression(10^1), expression(10^2),
                    expression(10^3), expression(10^4),
                    expression(10^5), expression(10^6),
                    expression(10^7), expression(10^8),
                    expression(10^9)) )
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
                       2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,
                       2e7,3e7,4e7,5e7,6e7,7e7,8e7,9e7) ), col = rgb(0,0,0,0.1) )
  abline( v = log10( c(1e1,2e1,3e1,4e1,5e1,6e1,7e1,8e1,9e1,
                       1e2,2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,
                       1e3,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,
                       1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,
                       1e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,
                       1e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,
                       1e7,2e7,3e7,4e7,5e7,6e7,7e7,8e7,9e7,
                       1e8,2e8,3e8,4e8,5e8,6e8,7e8,8e8,9e8,
                       1e9,2e9,3e9,4e9,5e9,6e9,7e9,8e9,9e9) ), col = rgb(0,0,0,0.1) )
  polygon( x = c( pyro_q$discharge, rev(pyro_q$discharge) ),
           y = c( pyro_q$min, rev(pyro_q$max) ),
           col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
  polygon( x = c( pyro_q$discharge, rev(pyro_q$discharge) ),
           y = c( pyro_q$lwr_iqr, rev(pyro_q$upr_iqr) ),
           col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
  # segments( x0 = 9.3,
  #           y0 = 0, y1 = pyro_q$upr_iqr[which(pyro_q$discharge==9.3)],
  #           lty = 2, col = 2 )
  # text( x = 9.3, y = pyro_q$upr_iqr[which(pyro_q$discharge==9.3)],
  #       col = 2, labels = "9.3 ug/L", pos = 4, srt = 90 )
  # segments( x0 = 0, x1 = 9.3,
  #           y0 = pyro_q$median[which(pyro_q$discharge==9.3)],
  #           lty = 2, col = 2 )
  # segments( x0 = 0, x1 = 9.3,
  #           y0 = pyro_q$lwr_iqr[which(pyro_q$discharge==9.3)],
  #           lty = 2, col = 2 )
  # segments( x0 = 0, x1 = 9.3,
  #           y0 = pyro_q$upr_iqr[which(pyro_q$discharge==9.3)],
  #           lty = 2, col = 2 )
  legend( 'bottomright', bty = 'n',
          legend = c("Median of maxima",
                     "IQR of maxima",
                     "Range of maxima (min/max)"), text.font = 2,
          text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )
  
  # 
  # #
  # #
  # #  There should be a single plot for annual mean(sd) discharge,
  # #  not separate plots for each subsegment !
  # #
  # #
  # 
  # # Plot monthly discharge (with 1sd and 2sd bands) by year, for reference
  # # initialize dataframe
  # q_year <- data.frame( year = 2012:2023,
  #                       mean = NA,
  #                       lwr_1sd = NA,
  #                       upr_1sd = NA,
  #                       lwr_2sd = NA,
  #                       upr_2sd = NA )
  # 
  # # populate rows
  # for( i in 1:nrow(q_year) ){
  #   # subset discharge data by year
  #   this <- pcdat$discharge[ which( year(pcdat$month)==q_year$year[i] ) ]
  #   # add output to q_year dataframe
  #   q_year[i,2:6] <- data.frame( mean = mean(this),
  #                                lwr_1sd = mean(this)-sd(this),
  #                                upr_1sd = mean(this)+sd(this),
  #                                lwr_2sd = mean(this)-2*sd(this),
  #                                upr_2sd = mean(this)+2*sd(this)
  #   )
  # }  # // end i loop
  # 
  # # Plot annual discharge values with 9.3 ug/L threshold
  # plot( mean ~ year, data = q_year, type = 'l', las = 1,
  #       main = paste0( subseg, " OTB\nMonthly discharge (non-zero Pyro months)" ),
  #       ylim = c(0,10), yaxt = 'n',
  #       ylab = "Discharge (m3/month)", xlab = "",
  #       lwd = 2, col = rgb(0,0.2,0.6,0.7))
  # axis( 2, at = 0:10, las = 1,
  #       labels = c( 0,expression(10^1), expression(10^2),
  #                   expression(10^3), expression(10^4),
  #                   expression(10^5), expression(10^6),
  #                   expression(10^7), expression(10^7),
  #                   expression(10^9), expression(10^10)) )
  # abline( h = seq(0,10,1), col = rgb(0,0,0,0.1) )
  # abline( h = log10( c(2e0,3e0,4e0,5e0,6e0,7e0,8e0,9e0,
  #                      2e1,3e1,4e1,5e1,6e1,7e1,8e1,9e1,
  #                      2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,
  #                      2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,
  #                      2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,
  #                      2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,
  #                      2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,
  #                      2e7,3e7,4e7,5e7,6e7,7e7,8e7,9e7,
  #                      2e8,3e8,4e8,5e8,6e8,7e8,8e8,9e8,
  #                      2e9,3e9,4e9,5e9,6e9,7e9,8e9,9e9,
  #                      2e10,3e10,4e10,5e10,6e10,7e10,8e10,9e10) ), col = rgb(0,0,0,0.1) )
  # abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  # polygon( x = c( q_year$year, rev(q_year$year) ),
  #          y = c( q_year$lwr_1sd, rev(q_year$upr_1sd) ),
  #          col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
  # polygon( x = c( q_year$year, rev(q_year$year) ),
  #          y = c( q_year$lwr_2sd, rev(q_year$upr_2sd) ),
  #          col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
  # # abline( h = 9.3, lty = 2, col = 2, lwd = 1.5 )
  # # text( x = 2022, y = 9.3,
  # #       col = 2, labels = "9.3 ug/L", pos = 3 )
  # points( mean ~ year, data = q_year )
  # legend( 'topright', bty = 'n',
  #         legend = c("mean","\u00B11sd","\u00B12sd"), text.font = 2,
  #         text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )
  
}  # end subseg loop

dev.off()
