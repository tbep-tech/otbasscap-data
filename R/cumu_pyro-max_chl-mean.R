rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )
load( "../data/Pyro.Rdata")

# Subset routine Pyro samples
pyro <- pyro[ which(pyro$routine==TRUE), ]

# Specify subsegment
subsegs <- c("NW","NE","CW","CE","SW","SE")

# Loop over sub-segments to generate plots
png( "../figs/cumu_pyro-max_chl-mean.png",
     width = 8, height = 10, units = 'in', res = 500 )
par( mfrow=c(3,2), mar=c(4,4,2.5,1) )
for( subseg in subsegs ){

  # Assemble data for analysis
    # chlorophyll data (EPC)
    epcwq3.sub <- epcwq3[ which( epcwq3$param=="Chla" &
                                 year(epcwq3$date) >= 2012 &
                                 epcwq3$subseg==subseg ), ]
    epcwq3.sub$month <- floor_date( epcwq3.sub$date, unit = 'month' ) 
    chldat <- epcwq3.sub |> group_by(month) |> dplyr::summarise( chl = mean(value) ) |> as.data.frame()
    # pyro data (FWC)
    pyro.sub <- pyro[ which( pyro$subsegment==subseg & pyro$yr>=2012 ), ]
    pyro.sub$month <- floor_date( pyro.sub$date, unit = 'month' ) 
    pyro.sub <- pyro.sub[ which(complete.cases(pyro.sub)), ]  # pyro==NA means zero cells/L
    pyro.sub$logval <- log10( pyro.sub$pyro )
    pyrodat <- pyro.sub |> group_by(month) |> dplyr::summarise( pyro = max(logval) ) |> as.data.frame()
    # join pyro and chl data by month
    pcdat <- inner_join( pyrodat, chldat, by = 'month' )
  
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
    plot( chl ~ median, data = pyro_chl, las = 1,
          type = 'l', lwd = 2, col = rgb(0,0.2,0.6,0.7),
          main = "", xlab = "", ylab = "",
          xlim = c(1,7), xaxt = 'n'
          )
    mtext( paste0(subseg, " sub-segment"),
           side = 3, adj = 0, line = 1 )
    # mtext( "Conditional distribution plot",
    #        side = 3, adj = 0, line = 0.1 )
    mtext( expression(italic(P.~bahamense)*" (cells/L)"),
           side = 1, line = 2.5, cex = 0.8 )
    mtext( "Chlorophyll a (ug/L)", side = 2, line = 2.5, cex = 0.8 )
    axis( 1, at = 1:7, las = 1,
          labels = c( expression(10^1), expression(10^2),
                      expression(10^3), expression(10^4),
                      expression(10^5), expression(10^6),
                      expression(10^7) ) )
    abline( v = log10( c(1e0,2e0,3e0,4e0,5e0,6e0,7e0,8e0,9e0,
                         1e1,2e1,3e1,4e1,5e1,6e1,7e1,8e1,9e1,
                         1e2,2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,
                         1e3,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,
                         1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,
                         1e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,
                         1e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,
                         1e7,2e7,3e7,4e7,5e7,6e7,7e7,8e7,9e7) ), col = rgb(0,0,0,0.1) )
    abline( h = min(pyro_chl$chl):max(pyro_chl$chl), col = rgb(0,0,0,0.1) )
    polygon( x = c( pyro_chl$min, rev(pyro_chl$max) ),
             y = c( pyro_chl$chl, rev(pyro_chl$chl) ),
             col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
    polygon( x = c( pyro_chl$lwr_iqr, rev(pyro_chl$upr_iqr) ),
             y = c( pyro_chl$chl, rev(pyro_chl$chl) ),
             col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
    segments( x0 = 0, x1 = pyro_chl$upr_iqr[which(pyro_chl$chl==9.3)],
              y0 = 9.3,
              lty = 3, col = rgb(1,0.1,0.1,0.6), lwd = 2 )
    text( x = 2, y = 9, col = rgb(1,0.1,0.1,0.8),
          labels = "9.3 ug/L", pos = 3, cex = 1.3 )
    segments( x0 = pyro_chl$median[which(pyro_chl$chl==9.3)],
              y0 = 0, y1 = 9.3,
              lty = 1, col = rgb(1,0.1,0.1,0.6), lwd = 2 )
    segments( x0 = pyro_chl$lwr_iqr[which(pyro_chl$chl==9.3)],
              y0 = 0, y1 = 9.3,
              lty = 3, col = rgb(1,0.1,0.1,0.6), lwd = 2 )
    segments( x0 = pyro_chl$upr_iqr[which(pyro_chl$chl==9.3)],
              y0 = 0, y1 = 9.3,
              lty = 3, col = rgb(1,0.1,0.1,0.6), lwd = 2 )
    legend( x = 0.5, y = 26, bty = 'n',
            legend = c("Median",
                       "IQR",
                       "Range"), text.font = 2, cex = 1.3,
            text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )
    

  # # Plot annual mean chlorophyll conc (with 1sd and 2sd bands) for reference
  #   # initialize dataframe
  #   chl_year <- data.frame( year = 2012:2023,
  #                           mean = NA,
  #                           lwr_1sd = NA,
  #                           upr_1sd = NA,
  #                           lwr_2sd = NA,
  #                           upr_2sd = NA )
  #   # populate rows
  #   for( i in 1:nrow(chl_year) ){
  #     # subset chl data by year
  #     this <- pcdat$chl[ which( year(pcdat$month)==chl_year$year[i] ) ]
  #     # add output to chl_year dataframe
  #     chl_year[i,2:6] <- data.frame( mean = mean(this),
  #                                    lwr_1sd = mean(this)-sd(this),
  #                                    upr_1sd = mean(this)+sd(this),
  #                                    lwr_2sd = mean(this)-2*sd(this),
  #                                    upr_2sd = mean(this)+2*sd(this)
  #                                    )
  #   }  # // end i loop
  # 
  #   # Plot annual mean chl values with 9.3 ug/L threshold
  #   plot( mean ~ year, data = chl_year, type = 'l', las = 1,
  #         ylim = c(0,50),
  #         main = "", ylab = "", xlab = "",
  #         cex.lab = 1.3,
  #         lwd = 2, col = rgb(0,0.2,0.6,0.7))
  #   # mtext( paste0(subseg, " sub-segment"),
  #   #        side = 3, adj = 0, line = 1.25 )
  #   mtext( "Annual mean chlorophyll-a",
  #          side = 3, adj = 0, line = 0.1 )
  #   mtext( "Chlorophyll a (ug/L)", side = 2, line = 2.5, cex = 0.8 )
  #   abline( h = seq(0,50,5), col = rgb(0,0,0,0.1) )
  #   abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  #   polygon( x = c( chl_year$year, rev(chl_year$year) ),
  #            y = c( chl_year$lwr_1sd, rev(chl_year$upr_1sd) ),
  #            col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
  #   polygon( x = c( chl_year$year, rev(chl_year$year) ),
  #            y = c( chl_year$lwr_2sd, rev(chl_year$upr_2sd) ),
  #            col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
  #   abline( h = 9.3, lty = 2, col = rgb(1,0.2,0.1,0.7), lwd = 2 )
  #   # text( x = 2022, y = 9.3,
  #   #       col = 2, labels = "9.3 ug/L", pos = 3 )
  #   points( mean ~ year, data = chl_year )
  #     legend( 'topright', bty = 'n',
  #             legend = c("Mean","\u00B11sd","\u00B12sd"), text.font = 2, cex = 1.3,
  #             text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.6), rgb(0,0.2,0.6,0.4) ) )

}  # end subseg loop

dev.off()
