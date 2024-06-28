rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data/Pyro.Rdata")

# Process pyro data
pyro$pyro[ which(is.na(pyro$pyro)) ] <- 0
pyro$logcells <- log10( pyro$pyro + 1 )
pyro <- pyro[ which( pyro$yr > 2011 ), ]

# Specify threshold for bloom events
threshold.low <- 1e4  # cells/L
threshold.med <- 1e5  
threshold.hi  <- 1e6


# OTB (aggregate) ---------------------------------------------------------

  # Plotting parameters
  png( "../figs/pyro_events_otb.png", height = 7, width = 10, units = 'in', res = 500 )
  layout( matrix(c(1,1,1, 2,3,4 ), 2, 3, byrow = TRUE) )
  par(mar=c(4,4,2,1))

  # A. Pyro sample data over time
  plot( logcells ~ date, data = pyro, col = rgb(0,0,0,0),
        ylim = c(0,7), las = 1,
        main = "(a) Pyrodinium abundance",
        ylab = "cells/L", xlab = "", yaxt = 'n', xaxt = 'n',
        cex.main = 1.3, cex.lab = 1.3 )
  axis( 1, at = as.Date( paste0( min(pyro$yr):(max(pyro$yr)+1),"-01-01" ) ),
        labels = min(pyro$yr):(max(pyro$yr)+1),
        cex.axis = 1.3 )
  axis( 2, at = 0:7, las = 1, cex.axis = 1.3,
        labels = c( 0, expression("10"^1), expression("10"^2),
                    expression("10"^3), expression("10"^4),
                    expression("10"^5), expression("10"^6), 
                    expression("10"^7) ) )
  abline( v = as.Date( paste0( min(pyro$yr):(max(pyro$yr)+1),"-01-01" ) ),
          col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  abline( h = log10(c(threshold.low,threshold.med,threshold.hi)),
          lty = 2, lwd = 1.5, col = rgb(0,0,0,0.6) )
  segments( x0 = pyro$date,
            y0 = rep( 0, nrow(pyro) ),
            y1 = pyro$logcells,
            col = rgb(0.9,0.1,0.2,0.8)
            )
  text( x = max(pyro$date),
        y = log10(c(threshold.low,threshold.med,threshold.hi)),
        labels = c("Low","Med","High"),
        pos = 3, col = rgb(0,0,0,0.8) )

  
  # B. Pyro statistics per year
  # Initiate dataframe
  pyro_stats <- data.frame( year = min(pyro$yr):max(pyro$yr),
                            min = NA,
                            lwr_iqr = NA,
                            median = NA,
                            upr_iqr = NA,
                            max = NA )
  # Populate rows
  for( i in 1:nrow(pyro_stats) ){
    this <- pyro$logcells[ which( pyro$yr==pyro_stats$year[i] &
                                  pyro$pyro > 0 ) ]
    pyro_stats$min[i] <- min(this)
    pyro_stats$lwr_iqr[i] <- quantile(this,0.25)
    pyro_stats$median[i] <- median(this)
    pyro_stats$upr_iqr[i] <- quantile(this,0.75)
    pyro_stats$max[i] <- max(this)
  }  # // end i loop
  # Plot
  plot( median ~ year, data = pyro_stats,
        col = rgb(0,0,0,0), type = 'l', ylim = c(0,7),
        main = "(b) Annual bloom statistics",
        ylab = "cells/L", xlab = "", yaxt = 'n',
        cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3
        )
  axis( 2, at = 0:7, las = 1, cex.axis = 1.3,
        labels = c( 0, expression("10"^1), expression("10"^2),
                    expression("10"^3), expression("10"^4),
                    expression("10"^5), expression("10"^6), 
                    expression("10"^7) ) )
  abline( v = pyro_stats$year, col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  polygon( x = c( pyro_stats$year, rev(pyro_stats$year) ),
           y = c( pyro_stats$min, rev(pyro_stats$max) ),
           border = rgb(0,0,0,0), col = rgb(1,0.1,0.2,0.2)
           )
  polygon( x = c( pyro_stats$year, rev(pyro_stats$year) ),
           y = c( pyro_stats$lwr_iqr, rev(pyro_stats$upr_iqr) ),
           border = rgb(0,0,0,0), col = rgb(1,0.1,0.2,0.3)
           )
  lines( median ~ year, data = pyro_stats, col = 2, lwd = 3 )
  points( median ~ year, data = pyro_stats, pch = 21,
          col = 2, bg = rgb(1,1,1,1), lwd = 2 )  
  legend( 'bottomleft', bty = 'n',
          legend = c("Min/max", "IQR", "Median"),
          text.col = c(rgb(1,0.1,0.2,0.4),rgb(1,0.1,0.2,0.7),2),
          text.font = 2 )
  legend( 'bottomright', bty = 'n',
          legend = "Null samples\nomitted\n",
          text.col = rgb(0,0,0,0.6), text.font = 2 )

  
  # C. Pyro bloom event start and end dates
  # Initialize dataframe
  events <- data.frame( year = min(pyro$yr):max(pyro$yr),
                        str_low = NA,
                        end_low = NA,
                        str_med = NA,
                        end_med = NA,
                        str_hi = NA,
                        end_hi = NA
                        )
  # Populate rows
  for( i in 1:nrow(events) ){
    event.dates.low <- pyro$doy[ which( pyro$yr==events$year[i] &
                                        pyro$pyro >= threshold.low ) ] |> sort()
    event.dates.med <- pyro$doy[ which( pyro$yr==events$year[i] &
                                        pyro$pyro >= threshold.med ) ] |> sort()
    event.dates.hi  <- pyro$doy[ which( pyro$yr==events$year[i] &
                                        pyro$pyro >= threshold.hi ) ] |> sort()
    events$str_low[i] <- event.dates.low |> min()
    events$end_low[i] <- event.dates.low |> max()
    events$str_med[i] <- event.dates.med |> min()
    events$end_med[i] <- event.dates.med |> max()
    events$str_hi[i]  <- event.dates.hi  |> min()
    events$end_hi[i]  <- event.dates.hi  |> max()
  }  # // end i loop
  # Manual corrections
  events[ which(events$year==2016), c("str_hi","end_hi") ] <- c(NA,NA)
  # Compute bloom durations by category
  events$duration_low <- events$end_low - events$str_low + 1
  events$duration_med <- events$end_med - events$str_med + 1
  events$duration_hi  <- events$end_hi - events$str_hi + 1
  # Manual corrections
  events$duration_hi[ which(events$year==2016) ] <- 0
  # Plot
  plot( str_med ~ year, data = events,
        type = 'l', lwd = 2, col = rgb(0,0,0,0),
        ylim = c(1,365), yaxt = 'n',
        main = "(c) Bloom event dates",
        ylab = "Day of year", xlab = "",
        cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3
        )
  axis( 2, at = seq(0,365,60), labels = seq(0,365,60),
        cex.axis = 1.3, las = 1 )
  abline( h = seq(0,365,60), col = rgb(0,0,0,0.1) )
  abline( v = events$year, col = rgb(0,0,0,0.1) )
  for( i in 1:nrow(events) ){
    polygon( x = c( rep( events$year[i]-0.12, 2), 
                    rep( events$year[i]+0.12, 2) ),
             y = c( events$str_low[i], events$end_low[i],
                    events$end_low[i], events$str_low[i]),
             col = rgb(0.3,0.1,0.2,0.3), border = rgb(1,1,1,0.4)
             )
    polygon( x = c( rep( events$year[i]-0.12, 2), 
                    rep( events$year[i]+0.12, 2) ),
             y = c( events$str_med[i], events$end_med[i],
                    events$end_med[i], events$str_med[i]),
             col = rgb(1,0.1,0.2,0.4), border = rgb(1,1,1,0.4)
             )
    if( events$duration_hi[i]>1 ){
      polygon( x = c( rep( events$year[i]-0.18, 2), 
                      rep( events$year[i]+0.18, 2) ),
               y = c( events$str_hi[i], events$end_hi[i],
                      events$end_hi[i], events$str_hi[i]),
               col = rgb(1,0.1,0.2,0.6), border = rgb(1,1,1,0.4)
               )
    } else if( events$duration_hi[i]==1 ) {
      segments( x0 = c(events$year[i]-0.18),
                x1 = c(events$year[i]+0.18),
                y0 = events$str_hi[i],
                col = rgb(1,0.1,0.2,0.6)
                )
    }
    
  }  # end i loop
  legend( 'bottomleft', bty = 'n',
          legend = c( expression(bold("Low (>10"^4*" cells/L)")),
                      expression(bold("Med (>10"^5*" cells/L)")),
                      expression(bold("High (>10"^6*" cells/L)")) ),
          text.col = c(rgb(0.3,0.1,0.2,0.5),rgb(1,0.1,0.2,0.6),
                       rgb(1,0.1,0.2,0.9)),
          text.font = 2 )
  
  
  # D. Bloom duration over time
  plot( duration_hi ~ year, data = events, las = 1,
        ylim = c(0,180), yaxt = 'n',
        type = 'l', col = rgb(0,0,0,0), lwd = 2,
        xlab = "", ylab = "Duration (days)",
        main = "(d) Bloom event duration",
        cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3
        )
  axis( 2, at = seq(0,180,30), labels = seq(0,180,30),
        cex.axis = 1.3, las = 1 )
  abline( h = seq(0,180,30), col = rgb(0,0,0,0.1) )
  abline( v = events$year, col = rgb(0,0,0,0.1) )
  polygon( x = c( events$year, rev(events$year)),
           y = c( events$duration_low, rep(0,nrow(events)) ),
           border = rgb(0,0,0,0), col = rgb(1,0.1,0.2,0.2)
           )
  polygon( x = c( events$year, rev(events$year)),
           y = c( events$duration_med, rep(0,nrow(events)) ),
           border = rgb(0,0,0,0), col = rgb(1,0.1,0.2,0.3)
           )
  polygon( x = c( events$year, rev(events$year)),
           y = c( events$duration_hi, rep(0,nrow(events)) ),
           border = rgb(0,0,0,0), col = rgb(1,0.1,0.2,0.3)
           )
  legend( 'topleft', bty = 'n',
          legend = c( expression(bold("Low (>10"^4*" cells/L)")),
                      expression(bold("Med (>10"^5*" cells/L)")),
                      expression(bold("High (>10"^6*" cells/L)")) ),
          text.col = c(rgb(1,0.1,0.2,0.4),rgb(1,0.1,0.2,0.6),
                       rgb(1,0.1,0.2,0.9)),
          text.font = 2 )
  
  dev.off()  # close png device
  
  
  # Additional plots
  #
  par( mfrow=c(2,2) )

  # E. Proportion of samples above each threshold, per event
    # Initialize dataframe
    pyro_prop <- data.frame( year = min(pyro$yr):max(pyro$yr),
                             prop_low = NA,
                             prop_med = NA,
                             prop_hi  = NA
    )
    # Populate rows
    for( i in 1:nrow(pyro_prop) ){
      this <- pyro$logcells[ which(pyro$yr==pyro_prop$year[i]) ]
      pyro_prop$prop_low[i] <- length(which(this>log10(threshold.low))) / length(this)
      pyro_prop$prop_med[i] <- length(which(this>log10(threshold.med))) / length(this)
      pyro_prop$prop_hi[i]  <- length(which(this>log10(threshold.hi)))  / length(this)
    }  # // end i loop
    # Plot
    plot( prop_low ~ year, data = pyro_prop, las = 1,
          col = rgb(0,0,0,0), ylim = c(0,0.40),
          main = "(e) Proportion of samples in each bloom category",
          xlab = "", ylab = "Proportion",
          cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3
    )
    abline( h = seq(0,1,0.05), col = rgb(0,0,0,0.1) )
    abline( v = pyro_prop$year, col = rgb(0,0,0,0.1) )
    
    polygon( x = c( pyro_prop$year, rev(pyro_prop$year)),
             y = c( pyro_prop$prop_low, rep(0,nrow(pyro_prop)) ),
             border = rgb(0,0,0,0), col = rgb(1,0.1,0.2,0.2)
    )
    polygon( x = c( pyro_prop$year, rev(pyro_prop$year)),
             y = c( pyro_prop$prop_med, rep(0,nrow(pyro_prop)) ),
             border = rgb(0,0,0,0), col = rgb(1,0.1,0.2,0.3)
    )
    polygon( x = c( pyro_prop$year, rev(pyro_prop$year)),
             y = c( pyro_prop$prop_hi, rep(0,nrow(pyro_prop)) ),
             border = rgb(0,0,0,0), col = rgb(1,0.1,0.2,0.3)
    )
    
    legend( 'topleft', bty = 'n',
            legend = c( expression(bold("Low (>10"^4*" cells/L)")),
                        expression(bold("Med (>10"^5*" cells/L)")),
                        expression(bold("High (>10"^6*" cells/L)")) ),
            text.col = c(rgb(1,0.1,0.2,0.4),rgb(1,0.1,0.2,0.7),
                         rgb(1,0.1,0.2,0.9)),
            text.font = 2 )
  
  # F. Proportion-duration relationships
    # Assemble data
    pyro_prop2 <- left_join( pyro_prop,
                             events[,c("year","duration_low",
                                       "duration_med","duration_hi")],
                             by = "year" )
    # Plot
    plot( prop_hi ~ duration_hi, data = pyro_prop2,
          xlim = c(0,150), ylim = c(0,0.25), las = 1,
          pch = 1, col = rgb(1,0.1,0.2,0.8), cex = 1.4, lwd = 1.5,
          main = "(f) Category proportions vs. bloom duration",
          xlab = "Duration (days)", ylab = "Proportion",
          cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3
    )
    abline( h = seq(0,1,0.05), col = rgb(0,0,0,0.1) )
    abline( v = seq(0,200,25), col = rgb(0,0,0,0.1) )
    lm_hi <- lm( prop_hi ~ duration_hi, data = pyro_prop2 )
    segments( x0 = min(pyro_prop2$duration_hi),
              x1 = max(pyro_prop2$duration_hi),
              y0 = lm_hi$coefficients[1] + lm_hi$coefficients[2]*min(pyro_prop2$duration_hi),
              y1 = lm_hi$coefficients[1] + lm_hi$coefficients[2]*max(pyro_prop2$duration_hi),
              lwd = 1.5, lty = 3
    )
    text( x = 80, y = 0.03,
          cex = 1, font = 2, col = rgb(0,0,0,0.6), adj = 0,
          labels = paste0( "High (>1M cells/L)",
                           "\nR2 = ", round( summary(lm_hi)$r.squared, 3),
                           "\nP < 0.001" )
    )
    points( prop_med ~ duration_med, data = pyro_prop2,
            pch = "+", col = rgb(1,0.1,0.2,0.6), cex = 1.4, lwd = 1.5 )
    lm_med <- lm( prop_med ~ duration_med, data = pyro_prop2 )
    segments( x0 = min(pyro_prop2$duration_med),
              x1 = max(pyro_prop2$duration_med),
              y0 = lm_med$coefficients[1] + lm_med$coefficients[2]*min(pyro_prop2$duration_med),
              y1 = lm_med$coefficients[1] + lm_med$coefficients[2]*max(pyro_prop2$duration_med),
              lwd = 1.5, lty = 3
    )
    text( x = 80, y = 0.17,
          cex = 1, font = 2, col = rgb(0,0,0,0.6), adj = 2, pos = 2,
          labels = paste0( "Med (>100k cells/L)",
                           "\nR2 = ", round( summary(lm_med)$r.squared, 3),
                           "\nP < 0.01" )
    )
  
  # G. Annual sample counts
    # Assemble sample counts (all samples and nonzero samples)
    pyro.N <- table( pyro$yr ) |> as.data.frame()
    colnames( pyro.N ) <- c("year","All samples")
    pyro.N$year <- pyro.N$year |> as.character() |> as.integer()
    pyro.N$`Nonzero samples` <- as.data.frame( table( pyro$yr[which(pyro$pyro>0)] ) )[,2]
    # Plot
    barplot( t(`row.names<-`(as.matrix(pyro.N[-1]), pyro.N$year)),
             legend = TRUE, beside = TRUE,
             args.legend = list( bty = 'n', x = 'top' ),
             main = "(g) Pyrodinium sample counts",
             ylab = "Count",
             cex.lab = 1.3, cex.main = 1.3, cex.axis = 1.3 )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
    
  # H. Proportion of non-zero samples
    # Compute proportion
    pyro.N$p <- pyro.N$`Nonzero samples` / pyro.N$`All samples`
    # Plot
    plot( p ~ year, data = pyro.N, las = 1,
          type = 'l', col = 2, lwd = 3,
          xlab = "", ylab = "Proportion",
          ylim = c( min(pyro.N$p)*0.85, max(pyro.N$p)*1.1 ),
          main = "(h) Proportion of non-zero samples",
          cex.lab = 1.3, cex.main = 1.3, cex.axis = 1.3 )
    points( p ~ year, data = pyro.N,
            pch = 21, cex = 1.2, lwd = 2, col = 2, bg = rgb(1,1,1,1) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )  
    
  

# OTB subsegments ---------------------------------------------------------

    # Get date of first observed 10^5 cell count each year
    med_str <- data.frame( yr = unique(pyro$yr),
                           str = as.Date(NA) )
    for( i in 1:nrow(med_str) ){
      med_str$str[i] <- pyro$date[ which( pyro$yr==med_str$yr[i] &
                                            pyro$logcells >= 5 ) ] |> min()
    }
    
    # Find subsegment(s) where 10^5 cell counts were first observed each year
    med_subseg <- list()
    for( i in 1:length(med_str$yr) ){
      med_subseg[[i]] <- pyro$subsegment[ which( pyro$yr == med_str$yr[i] &
                                            pyro$logcells >= 5 &
                                            floor_date(pyro$date,'week') == floor_date(med_str$str[i],'week') ) ] |> unique()
       names(med_subseg)[i] <- med_str$yr[i]
    }
    
  # Plotting parameters
  png( "../figs/pyro_events_otb-subsegments.png",
       height = 10, width = 8, units = 'in', res = 500 )
  layout( matrix( c(1,1,1,1,
                    2,2,2,2,
                    3,3,3,3,
                    4,4,4,4,
                    5,5,5,5,
                    6,6,6,6,
                    7), ncol=1 ) )
  par(mar=c(1,4,1.5,1))
  
  # Loop over segments to generate plots
  subsegs <- c("NW","NE","CW","CE","SW","SE")
  for( i in 1:length(subsegs) ){
    this <- pyro[ which(pyro$subsegment==subsegs[i]), ]
    plot( logcells ~ date, data = this, col = rgb(0,0,0,0),
          ylim = c(0,7), las = 1,
          main = "",
          ylab = "cells/L", xlab = "", yaxt = 'n', xaxt = 'n',
          cex.main = 1.3, cex.lab = 1.3 )
    for( j in 1:length(med_subseg) ){
      if( subsegs[i] %in% unlist(med_subseg[j]) ){
        jdates <- c( as.Date(paste0(names(med_subseg)[j],"-01-01")),
                     as.Date(paste0(names(med_subseg)[j],"-12-31")) )
        polygon( x = c( jdates[1], jdates[2], jdates[2], jdates[1] ),
                 y = c(-1,-1,8,8),
                 border = rgb(0,0,0,0), col = rgb(1,0,0,0.1)
                 )
      }
    }  # // end j loop
    
    quarters <- seq.Date( floor_date(min(pyro$date),'year'),
                          ceiling_date(max(pyro$date),'year'), 'quarter' )
    years <- seq.Date( floor_date(min(pyro$date),'year'),
                       ceiling_date(max(pyro$date),'year'), 'year' )
    if( i==length(subsegs) ){
      axis( 1, at = years,
            labels = min(pyro$yr):(max(pyro$yr)+1),
            cex.axis = 1.3 )
    } else {
      axis( 1, at = years,
            labels = rep("",length(years)),
            cex.axis = 1.3 )
    }
    mtext( paste0( "(", letters[i], ") ", subsegs[i], " sub-segment" ),
           side = 3, font = 2, adj = 0 )
    axis( 2, at = 0:7, las = 1, cex.axis = 1.3,
          labels = c( 0, expression("10"^1), expression("10"^2),
                      expression("10"^3), expression("10"^4),
                      expression("10"^5), expression("10"^6), 
                      expression("10"^7) ) )
    abline( v = years, col = rgb(0,0,0,0.1) )
    abline( v = quarters, col = rgb(0,0,0,0.1) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
    abline( h = log10(c(threshold.low,threshold.med,threshold.hi)),
            lty = 2, lwd = 1.5, col = rgb(0,0,0,0.6) )
    segments( x0 = this$date,
              y0 = rep( 0, nrow(this) ),
              y1 = this$logcells,
              col = rgb(0.9,0.1,0.2,0.8)
    )
    text( x = max(pyro$date),
          y = log10(c(threshold.low,threshold.med,threshold.hi))-0.2,
          labels = c("Low","Med","High"),
          pos = 3, col = rgb(0,0,0,0.8) )
    # segments( x0 = med_str$str, y0=-1, y1=0, col = rgb(0,0.3,0.8,0.7) )
    points( x = med_str$str, y = rep(-0.2,nrow(med_str)), pch = 2 )
    
  }  # // end plotting loop
  
  dev.off()
  
  
  