
rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

load('../data/loads.RData')  # monthly TN loads

TN <- loads[ which( year(loads$date)>=2000 &
                   loads$param=="TN load"),
             c('date','value') ]
TN$month <- month( TN$date )
TN_mo <- TN |> group_by(month) |> dplyr::summarise(mean=mean(value)) |> as.data.frame()
TN_mo$pct <- TN_mo$mean / sum(TN_mo$mean)
TN_mo$TMDL_portion <- TN_mo$pct * 486


# Plot average TN load per month associated with 486 tons (TMDL). This is what
# the entities can strive to get to collectively on a monthly basis.
png( "../figs/TMDL_transition.png", width = 8, height = 5, units = 'in', res = 600 )
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot( TMDL_portion ~ month, data = TN_mo, ylim = c(0,125),
      # main = "Average monthly load distribution by month",
      xlab = "", ylab = "TN load (tons/month)",
      type='l', xaxt = 'n', las = 1, lwd = 3, col = rgb(0.5,0.5,0.5,0.3) )
axis( 1, at = 1:12, labels = month(1:12,label=TRUE,abbr=TRUE) )
abline( h = seq(0,150,10), col = rgb(0,0,0,0.1) )
points( TMDL_portion ~ month, data = TN_mo, pch = 21,
        bg = rgb(1,1,1,1), col = rgb(0.5,0.5,0.5,1), lwd = 2, cex = 1.2 )
text( x = TN_mo$month, y = TN_mo$TMDL_portion,
      pos = 3, col = rgb(0.2,0.5,0.9,1), font = 2,
      labels = paste0(round(100*TN_mo$pct,1),"%") )
dev.off()

# Plot the average TN load percentages by month. These percentages can be
# used as a guide by each entity to plan their monthly reductions.
par(mfrow=c(1,1))
plot( pct ~ month, data = TN_mo,
      main = "Average monthly percentage of total annual load (2000-2021)",
      xlab = "", ylab = "Mean percentage",
      type='l', yaxt = 'n', xaxt = 'n', lwd = 3 )
axis( 1, at = 1:12, labels = month(1:12,label=TRUE,abbr=TRUE) )
axis( 2, at = seq(0,0.40,0.05), las = 1,
      labels = paste0(seq(0,0.40,0.05)*100,"%") )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
points( pct ~ month, data = TN_mo, pch = 21,
        bg = rgb(1,1,1,1), lwd = 2, cex = 1.2 )




