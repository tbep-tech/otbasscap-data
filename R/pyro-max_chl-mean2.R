rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )
load( "../data/Pyro.Rdata" )
load( "../data-clean/epcphyto.Rdata" )

# Subset routine Pyro samples
pyro <- pyro[ which(pyro$routine==TRUE), ]

# Specify subsegment
subsegs <- c("NW","NE","CW","CE","SW","SE")

# Loop over sub-segments to generate plots
png( "../figs/pyro-max_chl-mean2.png", width = 7, height = 9, units = 'in', res = 500 )
par( mfrow=c(3,2), mar=c(4,4,3,2) )
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
    # phytoplankton data (EPC)
    phyto.sub <- phyto[ which( year(phyto$date) >= 2012 &
                               year(phyto$date) <= 2023 &
                               phyto$subsegment==subseg ), ]
    # phyto.sub <- phyto.sub[ which( phyto.sub$family=="Bacillariaceae" |
    #                                phyto.sub$name=="Pseudo-nitzschia sp." |
    #                                phyto.sub$name=="Tripos hircus" ), ]
    phyto.sub <- phyto.sub[ which( phyto.sub$family=="Bacillariaceae" ), ]
    phyto.sub$month <- floor_date( phyto.sub$date, unit = 'month' )
    phytodat <- phyto.sub |> group_by(month) |>
                  dplyr::summarise( cells = max(cellcount) ) |> as.data.frame()
    phytodat$phyto <- log10( phytodat$cells + 1 )
    phytodat <- phytodat[, c('month','phyto') ]
    # join pyro and chl data by month
    pcdat <- inner_join( pyrodat, chldat, by = 'month' )
    pcdat <- inner_join( pcdat, phytodat, by = 'month' )
    # pcdat$phyto <- pcdat$phyto + pcdat$pyro
  
  # Define function to summarize pyro or phyto distribution
    distn <- function( x, max_chl, col ){
      # subset pyro data associated with chl values at or below max_chl
      this <- x[ which( x$chl <= max_chl ), col ]
      # assemble output statistics
      out <- data.frame( median = median(this),
                         min = min(this),
                         lwr_iqr = quantile(this, 0.25),
                         upr_iqr = quantile(this, 0.75),
                         max = max(this)
                         )
      return( out )
    }  # // end distn()
  
  
  # Calculate phyto distn statistics
    # initiate dataframe
    phyto_chl <- data.frame( chl = seq(4,25,0.1),
                            median = NA,
                            min = NA,
                            lwr_iqr = NA,
                            upr_iqr = NA,
                            max = NA
                            )
    # populate rows
    for( i in 1:nrow(phyto_chl) ){
      phyto_chl[i,2:6] <- distn( pcdat, max_chl = phyto_chl[i,1], col = 'phyto' )
    }  # // end i loop
    
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
      pyro_chl[i,2:6] <- distn( pcdat, max_chl = pyro_chl[i,1], col = 'pyro' )
    }  # // end i loop
  
  # Plot pyro distn as a function of chl_max
    plot( median ~ chl, data = pyro_chl, las = 1,
          type = 'l', lwd = 2, col = rgb(0,0.2,0.6,0.7),
          main = paste0( subseg, " OTB\nMax Pyro distn ~ Mean Chl-a" ),
          ylim = c(1,7), yaxt = 'n',
          ylab = "cells/L", xlab = ""
          )
    mtext( "Chlorophyll a (ug/L)", side = 1, line = 2, cex = 0.7 )
    axis( 2, at = 1:7, las = 1,
          labels = c( expression(10^1), expression(10^2),
                      expression(10^3), expression(10^4),
                      expression(10^5), expression(10^6),
                      expression(10^7) ) )
    abline( h = 1:7, col = rgb(0,0,0,0.1) )
    abline( v = min(pyro_chl$chl):max(pyro_chl$chl), col = rgb(0,0,0,0.1) )
    polygon( x = c( pyro_chl$chl, rev(pyro_chl$chl) ),
             y = c( pyro_chl$min, rev(pyro_chl$max) ),
             col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
    polygon( x = c( pyro_chl$chl, rev(pyro_chl$chl) ),
             y = c( pyro_chl$lwr_iqr, rev(pyro_chl$upr_iqr) ),
             col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
    segments( x0 = 9.3,
              y0 = 0, y1 = par("usr")[4],
              lty = 2, col = 2 )
    text( x = 9.3, y = 2,
          col = 2, labels = "9.3 ug/L", pos = 1, srt = 90 )
    # segments( x0 = 0, x1 = 9.3,
    #           y0 = pyro_chl$median[which(pyro_chl$chl==9.3)],
    #           lty = 2, col = 2 )
    # segments( x0 = 0, x1 = 9.3,
    #           y0 = pyro_chl$lwr_iqr[which(pyro_chl$chl==9.3)],
    #           lty = 2, col = 2 )
    # segments( x0 = 0, x1 = 9.3,
    #           y0 = pyro_chl$upr_iqr[which(pyro_chl$chl==9.3)],
    #           lty = 2, col = 2 )
    # legend( 'bottomleft', bty = 'n',
    #         legend = c("Median of maxima",
    #                    "IQR of maxima",
    #                    "Range of maxima (min/max)"), text.font = 2,
    #         text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )
  
  # Plot median phyto~chl for reference
    lines( median ~ chl, data = phyto_chl, lwd = 2, col = rgb(0,0,0,0.5) )
    legend( 'bottomright', bty = 'n', lwd = 2,
            legend = c("Pyrodinium",
                       "Bacillariaceae"),
            col = c(rgb(0,0.2,0.6,0.9), rgb(0,0,0,0.7) ) )
  
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
  #     this <- epcwq3.sub$value[ which( year(epcwq3.sub$month)==chl_year$year[i] ) ]
  #     # add output to chl_year dataframe
  #     chl_year[i,2:6] <- data.frame( mean = mean(this),
  #                                    lwr_1sd = mean(this)-sd(this),
  #                                    upr_1sd = mean(this)+sd(this),
  #                                    lwr_2sd = mean(this)-2*sd(this),
  #                                    upr_2sd = mean(this)+2*sd(this)
  #                                    )
  #   }  # // end i loop

    # # Plot annual mean chl values with 9.3 ug/L threshold
    # plot( mean ~ year, data = chl_year, type = 'l', las = 1,
    #       main = paste0( subseg, " OTB\nAnnual mean chlorophyll-a concentration" ),
    #       ylim = c(0,50),
    #       ylab = "Chlorophyll a (ug/L)", xlab = "",
    #       lwd = 2, col = rgb(0,0.2,0.6,0.7))
    # abline( h = seq(0,50,5), col = rgb(0,0,0,0.1) )
    # abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
    # polygon( x = c( chl_year$year, rev(chl_year$year) ),
    #          y = c( chl_year$lwr_1sd, rev(chl_year$upr_1sd) ),
    #          col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
    # polygon( x = c( chl_year$year, rev(chl_year$year) ),
    #          y = c( chl_year$lwr_2sd, rev(chl_year$upr_2sd) ),
    #          col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
    # abline( h = 9.3, lty = 2, col = 2, lwd = 1.5 )
    # text( x = 2022, y = 9.3,
    #       col = 2, labels = "9.3 ug/L", pos = 3 )
    # points( mean ~ year, data = chl_year )
    # legend( 'topright', bty = 'n',
    #         legend = c("mean","\u00B11sd","\u00B12sd"), text.font = 2,
    #         text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )
    
    # # Plot phyto ~ pyro
    # # Assemble dataframe
    # phyto_pyro <- data.frame( chl = pyro_chl$chl,
    #                           pyro = pyro_chl$median )
    # phyto_pyro <- left_join( phyto_pyro, phyto_chl[,c("chl","median")], by = "chl" )
    # colnames(phyto_pyro)[3] <- "phyto"
    # phyto_pyro$col <- paste0( "rgb(1,0.2,0.1,",seq(0.15,0.80,length.out=nrow(phyto_pyro)),")" )
    # plot( phyto ~ pyro, data = phyto_pyro, pch = 16,
    #       col = apply( matrix(phyto_pyro$col,ncol=1), 1, function(x) eval(parse(text=x)) ),
    #       xlim = c(2,6), ylim = c(6,10)
    #       )
    # abline( v=1:12,col=rgb(0,0,0,0.1) ); abline(h=1:12,col=rgb(0,0,0,0.1))

}  # end subseg loop

dev.off()
