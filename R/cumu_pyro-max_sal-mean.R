rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )
load( "../data/Pyro.Rdata")

# Specify subsegment
subsegs <- c("NW","NE","CW","CE","SW","SE")

# Loop over sub-segments to generate plots
png( "../figs/cumu_pyro-max_sal-mean.png", width = 14, height = 9, units = 'in', res = 500 )
par( mfrow=c(3,4), mar=c(4,4,3,2) )
for( subseg in subsegs ){

  # Assemble data for analysis
    # Salinity data (EPC)
    epcwq3.sub <- epcwq3[ which( epcwq3$param=="Sal_top" &
                                 year(epcwq3$date) >= 2012 &
                                 epcwq3$subseg==subseg ), ]
    epcwq3.sub$month <- floor_date( epcwq3.sub$date, unit = 'month' ) 
    saldat <- epcwq3.sub |> group_by(month) |> summarise( sal = mean(value) ) |> as.data.frame()
    # pyro data (FWC)
    pyro.sub <- pyro[ which( pyro$subsegment==subseg & pyro$yr>=2012 ), ]
    pyro.sub$month <- floor_date( pyro.sub$date, unit = 'month' ) 
    pyro.sub <- pyro.sub[ which(complete.cases(pyro.sub)), ]  # pyro==NA means zero cells/L
    pyro.sub$logval <- log10( pyro.sub$pyro )
    pyrodat <- pyro.sub |> group_by(month) |> summarise( pyro = max(logval) ) |> as.data.frame()
    # join pyro and sal data by month
    wqdat <- inner_join( pyrodat, saldat, by = 'month' )
  
  # Define function to summarize pyro distribution
    distn <- function( x, max_sal ){
      # subset pyro data associated with sal values at or below max_sal
      this <- x$pyro[ which( x$sal <= max_sal ) ]
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
    pyro_sal <- data.frame( sal = seq(10,35,1),
                            median = NA,
                            min = NA,
                            lwr_iqr = NA,
                            upr_iqr = NA,
                            max = NA
                            )
    # populate rows
    for( i in 1:nrow(pyro_sal) ){
     pyro_sal[i,2:6] <- distn( wqdat, max_sal = pyro_sal[i,1] )
    }  # // end i loop
  
  # Plot pyro distn as a function of sal_max
    plot( median ~ sal, data = pyro_sal, las = 1,
          type = 'l', lwd = 2, col = rgb(0,0.2,0.6,0.7),
          main = paste0( subseg, " OTB\nMax Pyro distn ~ Mean Salinity" ),
          ylim = c(1,7), yaxt = 'n',
          ylab = "Pyro (cells/L)", xlab = ""
          )
    mtext( "Salinity (PSU)", side = 1, line = 2, cex = 0.7 )
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
    abline( v = min(pyro_sal$sal):max(pyro_sal$sal), col = rgb(0,0,0,0.1) )
    polygon( x = c( pyro_sal$sal, rev(pyro_sal$sal) ),
             y = c( pyro_sal$min, rev(pyro_sal$max) ),
             col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
    polygon( x = c( pyro_sal$sal, rev(pyro_sal$sal) ),
             y = c( pyro_sal$lwr_iqr, rev(pyro_sal$upr_iqr) ),
             col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
    # segments( x0 = 9.3,
    #           y0 = 0, y1 = pyro_sal$upr_iqr[which(pyro_sal$sal==9.3)],
    #           lty = 2, col = 2 )
    # text( x = 9.3, y = pyro_sal$upr_iqr[which(pyro_sal$sal==9.3)],
    #       col = 2, labels = "9.3 ug/L", pos = 4, srt = 90 )
    # segments( x0 = 0, x1 = 9.3,
    #           y0 = pyro_sal$median[which(pyro_sal$sal==9.3)],
    #           lty = 2, col = 2 )
    # segments( x0 = 0, x1 = 9.3,
    #           y0 = pyro_sal$lwr_iqr[which(pyro_sal$sal==9.3)],
    #           lty = 2, col = 2 )
    # segments( x0 = 0, x1 = 9.3,
    #           y0 = pyro_sal$upr_iqr[which(pyro_sal$sal==9.3)],
    #           lty = 2, col = 2 )
    legend( 'bottomright', bty = 'n',
            legend = c("Median of maxima",
                       "IQR of maxima",
                       "Range of maxima (min/max)"), text.font = 2,
            text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )
  
  
  # Plot annual mean salinity (with 1sd and 2sd bands) for reference
    # initialize dataframe
    sal_year <- data.frame( year = 2012:2023,
                            mean = NA,
                            lwr_1sd = NA,
                            upr_1sd = NA,
                            lwr_2sd = NA,
                            upr_2sd = NA )
    # populate rows
    for( i in 1:nrow(sal_year) ){
      # subset sal data by year
      this <- wqdat$sal[ which( year(wqdat$month)==sal_year$year[i] ) ]
      # add output to sal_year dataframe
      sal_year[i,2:6] <- data.frame( mean = mean(this),
                                     lwr_1sd = mean(this)-sd(this),
                                     upr_1sd = mean(this)+sd(this),
                                     lwr_2sd = mean(this)-2*sd(this),
                                     upr_2sd = mean(this)+2*sd(this)
                                     )
    }  # // end i loop
  
    # Plot annual mean sal values
    plot( mean ~ year, data = sal_year, type = 'l', las = 1,
          main = paste0( subseg, " OTB\nAnnual mean salinity (non-zero Pyro months)" ),
          ylim = c(0,40),
          ylab = "Salinity (PSU)", xlab = "",
          lwd = 2, col = rgb(0,0.2,0.6,0.7))
    abline( h = seq(0,50,5), col = rgb(0,0,0,0.1) )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
    polygon( x = c( sal_year$year, rev(sal_year$year) ),
             y = c( sal_year$lwr_1sd, rev(sal_year$upr_1sd) ),
             col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
    polygon( x = c( sal_year$year, rev(sal_year$year) ),
             y = c( sal_year$lwr_2sd, rev(sal_year$upr_2sd) ),
             col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
    # abline( h = 9.3, lty = 2, col = 2, lwd = 1.5 )
    # text( x = 2022, y = 9.3,
    #       col = 2, labels = "9.3 ug/L", pos = 3 )
    points( mean ~ year, data = sal_year )
    legend( 'topright', bty = 'n',
            legend = c("mean","\u00B11sd","\u00B12sd"), text.font = 2,
            text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )

}  # end subseg loop

dev.off()
