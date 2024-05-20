rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# load data
load(file = here::here('data-clean/loads.RData'))
pyro <- read.csv(file = here::here('data-raw/fwcpyro20112023.csv'))

# Pyrodinium
#
  pyro$date <- pyro$date |> as.Date(format='%Y-%m-%d')
  pyro$pyro[ which(is.na(pyro$pyro)) ] <- 0
  pyro$logcells <- log10( pyro$pyro + 1 )
  # Annual maximum Pyro abundance
  pyro.max <- pyro |> group_by(yr) |>
                dplyr::summarise( max = max(logcells,na.rm=TRUE) ) |>
                as.data.frame()
  pyro.max <- pyro.max[-1,]  # remove the 2011 record because data are sparse
  # Pyro bloom event duration
  threshold <- 1000
  pyro.dura <- pyro |> group_by(yr) |> 
                 dplyr::summarise( duration = as.integer(diff(range(date[which(pyro>threshold)]))) ) |> 
                 as.data.frame()
  pyro.dura <- pyro.dura[-1,]  # remove the 2011 record because data are sparse
    
  # Plots
  png("figs/Pyro_time.png", res = 400, units = 'in', width = 7, height = 10 )
  par(mfrow=c(3,1))
    # Pyro abundance (log scale)
    plot( logcells ~ date, data = pyro, col = rgb(0,0,0,0),
          ylim = c(0,7), las = 1,
          main = "Pyrodinium abundance",
          ylab = "cells/L", xlab = "", yaxt = 'n', xaxt = 'n',
          cex.main = 1.4, cex.lab = 1.3 )
    axis( 1, at = seq.Date(as.Date('2011-01-01'),
                           as.Date('2024-01-01'),'year'),
             labels = 2011:2024, cex.axis = 1.3 )
    axis( 2, at = 0:7, las = 1, cex.axis = 1.3,
          labels = c( 0, expression("10"^1), expression("10"^2),
                      expression("10"^3), expression("10"^4),
                      expression("10"^5), expression("10"^6), 
                      expression("10"^7) ) )
    abline( v = seq.Date(as.Date('2011-01-01'),
                         as.Date('2024-01-01'),'year'),
            col = rgb(0,0,0,0.2) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
    segments( x0 = pyro$date,
              y0 = rep( 0, nrow(pyro) ),
              y1 = pyro$logcells,
              col = rgb(0.7,0.3,0,1)
              )
    # Annual maximum Pyro abundance (log scale)
    plot( max ~ yr, data = pyro.max, type = 'l',
          lwd = 3, col = rgb(0.7,0.3,0,0.7), ylim = c(5,7),
          main = "Annual maximum Pyrodinium abundance",
          ylab = "log10( cells/L )", xlab = "", xaxt = 'n',
          cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.3, las = 1
          )
    axis( 1, at = 2011:2024,
          labels = 2011:2024, cex.axis = 1.3 )
    abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
    points( max ~ yr, data = pyro.max, pch = 21, lwd = 2.5, cex = 1.5,
            col = rgb(0.7,0.3,0,1), bg = rgb(1,1,1,1) )
    # Pyro bloom event duration abundance 
    plot( duration ~ yr, data = pyro.dura, type = 'l',
          lwd = 3, col = rgb(0.7,0.3,0,0.7), ylim = c(0,250),
          main = paste0("Pyro bloom event duration (min ",threshold," cells/L)"),
          ylab = "Duration (days)", xlab = "", xaxt = 'n',
          cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.3, las = 1
    )
    axis( 1, at = 2011:2024,
          labels = 2011:2024, cex.axis = 1.3 )
    abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
    points( duration ~ yr, data = pyro.dura, pch = 21, lwd = 2.5, cex = 1.5,
            col = rgb(0.7,0.3,0,1), bg = rgb(1,1,1,1) )

    dev.off()

# TN loads
#
  # Annual TN load to OTB
  loads$yr <- year( loads$date )
  load.annual <- loads |> group_by(yr) |> dplyr::summarise(load = sum(value)) |> 
                 as.data.frame()
  # "Wet season" TN load to OTB
  load.season <- loads[ which( month(loads$date)>3 & month(loads$date)<9 ), ] |>
                   group_by(yr) |>
                   dplyr::summarise(load = sum(value)) |> 
                   as.data.frame()
  
# Linear models
#
  png("figs/Pyro_TN.png", res = 400, units = 'in', width = 9, height = 9 )
  par(mfrow=c(2,2))
  # Annual max Pyro vs. annual TN load
  df1 <- inner_join( load.annual, pyro.max, by = 'yr' )
  mod1 <- lm( max ~ load, data = df1 )
  summary( mod1 )
  # Scatterplot + linear model
  plot( max ~ load, data = df1, pch = 21,
        main = "Annual max Pyro vs. annual TN load",
        xlab = "TN load, tons",
        ylab = "Maximum Pyrodinium abundance, log10(cells/L)",
        xlim = c(300,1000), ylim = c(5.5,7),
        cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.3,
        col = rgb(0.9,0.3,0,1), lwd = 2.5, cex = 1.5 )
  mtext( side = 3,
         text = paste0( "R2 = ",round(summary(mod1)$r.squared,3),
                        "; P = ",round(summary(mod1)$coefficients[2,4],3) ))
  abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
  segments( x0 = min(df1$load), x1 = max(df1$load),
            y0 = mod1$coefficients[1] + min(df1$load)*mod1$coefficients[2],
            y1 = mod1$coefficients[1] + max(df1$load)*mod1$coefficients[2],
            lty = 2, lwd = 3, col = rgb(0,0,0,0.8) )
  
  # Annual max Pyro vs. wet-season TN load
  df2 <- inner_join( load.season, pyro.max, by = 'yr' )
  mod2 <- lm( max ~ load, data = df2 )
  summary( mod2 )
  # Scatterplot + linear model
  plot( max ~ load, data = df2, pch = 21,
        main = "Annual max Pyro vs. wet-season TN load",
        xlab = "TN load, tons",
        ylab = "Maximum Pyrodinium abundance, log10(cells/L)",
        xlim = c(100,700), ylim = c(5.5,7),
        cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.3,
        col = rgb(0.9,0.3,0,1), lwd = 2.5, cex = 1.5 )
  mtext( side = 3,
         text = paste0( "R2 = ",round(summary(mod2)$r.squared,3),
                        "; P = ",round(summary(mod2)$coefficients[2,4],3) ))
  abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
  segments( x0 = min(df2$load), x1 = max(df2$load),
            y0 = mod2$coefficients[1] + min(df2$load)*mod2$coefficients[2],
            y1 = mod2$coefficients[1] + max(df2$load)*mod2$coefficients[2],
            lty = 2, lwd = 3, col = rgb(0,0,0,0.8) )

  # Annual max Pyro vs. annual TN load
  df3 <- inner_join( load.annual, pyro.dura, by = 'yr' )
  mod3 <- lm( duration ~ load, data = df3 )
  summary( mod3 )
  # Scatterplot + linear model
  plot( duration ~ load, data = df3, pch = 21,
        main = "Pyro bloom duration vs. annual TN load",
        xlab = "TN load, tons",
        ylab = "Pyrodinium bloom duration, days",
        xlim = c(300,1000), ylim = c(50,250),
        cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.3,
        col = rgb(0.9,0.3,0,1), lwd = 2.5, cex = 1.5 )
  mtext( side = 3,
         text = paste0( "R2 = ",round(summary(mod3)$r.squared,3),
                        "; P = ",round(summary(mod3)$coefficients[2,4],3) ))
  abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
  segments( x0 = min(df3$load), x1 = max(df3$load),
            y0 = mod3$coefficients[1] + min(df3$load)*mod3$coefficients[2],
            y1 = mod3$coefficients[1] + max(df3$load)*mod3$coefficients[2],
            lty = 2, lwd = 3, col = rgb(0,0,0,0.8) )
  
  # Annual max Pyro vs. wet-season TN load
  df4 <- inner_join( load.season, pyro.dura, by = 'yr' )
  mod4 <- lm( duration ~ load, data = df4 )
  summary( mod4 )
  # Scatterplot + linear model
  plot( duration ~ load, data = df4, pch = 21,
        main = "Pyro bloom duration vs. wet-season TN load",
        xlab = "TN load, tons",
        ylab = "Pyrodinium bloom duration, days",
        xlim = c(100,700), ylim = c(50,250),
        cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.3,
        col = rgb(0.9,0.3,0,1), lwd = 2.5, cex = 1.5 )
  mtext( side = 3,
         text = paste0( "R2 = ",round(summary(mod4)$r.squared,3),
                        "; P = ",round(summary(mod4)$coefficients[2,4],3) ))
  abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
  segments( x0 = min(df4$load), x1 = max(df4$load),
            y0 = mod4$coefficients[1] + min(df4$load)*mod4$coefficients[2],
            y1 = mod4$coefficients[1] + max(df4$load)*mod4$coefficients[2],
            lty = 2, lwd = 3, col = rgb(0,0,0,0.8) )
  
  dev.off()
  