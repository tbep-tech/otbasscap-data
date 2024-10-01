
rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

load('../data/loads.RData')  # monthly TN loads

# Compute mean & sd TN load by calendar month
TN <- loads[ which( year(loads$date)>=2000 &
                   loads$param=="TN load"),
             c('date','value') ]
TN$month <- month( TN$date )
TN_mo <- TN |> group_by(month) |> dplyr::summarise(mean=mean(value)) |> as.data.frame()
# TN_SE <- TN |> group_by(month) |> dplyr::summarise(SE=sd(value)/sqrt(length(unique(year(TN$date)))) ) |>
#            as.data.frame()

# TN_mo <- inner_join( TN_mo, TN_SE, by = "month" )
TN_mo$pct <- TN_mo$mean / sum(TN_mo$mean)

# Compute by-month TN loading consistent with absolute TMDL attainment
TN_mo$TMDL <- TN_mo$pct * 486

# Compute monthly cumulative loads
TN_mo$mean_cdf <- cumsum( TN_mo$mean )
TN_mo$TMDL_cdf <- cumsum( TN_mo$TMDL )

# Plot by-month TN load averages for comparison to absolute TMDL attainment
png( "../figs/loads_mo.png", width = 8, height = 8, units = 'in', res = 600 )
par(mfrow=c(2,1),mar=c(4,4,1,1))

  # TN load by month
  plot( TMDL ~ month, data = TN_mo, ylim = c(0,320),
        xlab = "", ylab = "TN load (tons/month)",
        type='l', xaxt = 'n', las = 1, lwd = 3, col = rgb(0,0,0,0) )
  mtext( "(a) TN load by month",
         side = 3, adj = 0, cex = 1.2 )
  polygon( x = c(5.5,9.5,9.5,5.5), y = c(-100,-100,500,500),
           border = rgb(0,0,0,0), col = rgb(0,0.5,0.8,0.1) )
  years <- TN$date |> year() |> unique() |> sort()
  for( i in 1:length(years) ){
    this.yr <- data.frame( month = TN$month[which(year(TN$date)==years[i]) ],
                            value = TN$value[which(year(TN$date)==years[i]) ] )
    lines( value ~ month, data = this.yr,
           lwd = 1, col = rgb(0,0.5,0.8,0.3) )
  }
  axis( 1, at = 1:12, labels = month(1:12,label=TRUE,abbr=TRUE) )
  abline( h = seq(0,600,25), col = rgb(0,0,0,0.1) )
  abline( h = seq(0,600,100), col = rgb(0,0,0,0.2) )
  lines( mean ~ month, data = TN_mo, lwd = 3, col = rgb(0,0.5,1,0.9) )
  points( mean ~ month, data = TN_mo, pch = 21,
          bg = rgb(1,1,1,1), col = rgb(0,0.5,1,1), lwd = 2, cex = 1.2 )
  lines( TMDL ~ month, data = TN_mo, lwd = 3, col = rgb(1,0.2,0.3,0.9) )
  points( TMDL ~ month, data = TN_mo, pch = 21,
          bg = rgb(1,1,1,1), col =  rgb(1,0.2,0.3,1), lwd = 2, cex = 1.2 )
  text( x = TN_mo$month, y = -20,
        pos = 3, col = rgb(0.3,0.3,0.4,1), font = 2, cex = 0.9,
        labels = paste0(format(round(100*TN_mo$pct,1),nsmall=1),"%") )
  legend( 'topleft', lwd = c(1,3,3), cex = 0.9,
          col = c(rgb(0,0.5,0.8,0.3),
                  rgb(0,0.5,1,0.9),
                  rgb(1,0.2,0.3,0.8)
          ),
          legend = c(paste0("TN load by month (2000","\U2013","2021)"),
                     paste0("Average TN load by month (2000","\U2013","2021)"),
                     'TMDL by month (total: 486 tons/year)'
          )
  )
  
  # Cumulative TN load by month
  plot( mean_cdf ~ month, data = TN_mo, ylim = c(0,900),
        xlab = "", ylab = "TN load (tons)",
        type='l', xaxt = 'n', las = 1, lwd = 3, col = rgb(0,0,0,0) )
  mtext( "(b) Cumulative TN load by month",
         side = 3, adj = 0, cex = 1.2 )
  polygon( x = c(5.5,9.5,9.5,5.5), y = c(-100,-100,1000,1000),
           border = rgb(0,0,0,0), col = rgb(0,0.5,0.8,0.1) )
  for( i in 1:length(years) ){
    this.cdf <- data.frame( month = TN$month[which(year(TN$date)==years[i]) ],
                            value = TN$value[which(year(TN$date)==years[i]) ] |> cumsum()
                           )
    lines( value ~ month, data = this.cdf,
           lwd = 1, col = rgb(0,0.5,0.8,0.3) )
  }
  axis( 1, at = 1:12, labels = month(1:12,label=TRUE,abbr=TRUE) )
  abline( h = seq(0,1000,50), col = rgb(0,0,0,0.1) )
  abline( h = seq(0,1000,100), col = rgb(0,0,0,0.2) )
  lines( mean_cdf ~ month, data = TN_mo, lwd = 3, col = rgb(0,0.5,1,0.9) )
  points( mean_cdf ~ month, data = TN_mo, pch = 21,
          bg = rgb(1,1,1,1), col = rgb(0,0.5,1,1), lwd = 2, cex = 1 )
  lines( TMDL_cdf ~ month, data = TN_mo, col = rgb(1,0.2,0.3,0.9), lwd = 3 )
  points( TMDL_cdf ~ month, data = TN_mo, col = rgb(1,0.2,0.3,1),
          pch = 21, bg = rgb(1,1,1,1), lwd = 2, cex = 1 )

dev.off()
