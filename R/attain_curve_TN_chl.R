rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# Load data
load('../data/loads.RData')  # monthly TN loads
load('../data-clean/epcwq_clean.RData')  # EPC wq data
load("../data-raw/totanndat.RData")

# Assemble annual chlorophyll attainment data
chldat <- epcwq3[ which( epcwq3$param=="Chla" &
                         year(epcwq3$date) >= 2000 ),
                  c('date','value') ]
chldat$month <- month( chldat$date )
chl_mo <- chldat |> group_by(date) |> dplyr::summarise(value=mean(value)) |> 
  as.data.frame()
chl_mo$year <- year( chl_mo$date )
chl_yr <- chl_mo |> group_by(year) |> dplyr::summarise(chl=mean(value)) |> 
  as.data.frame()
chl_yr$chl_attain <- ifelse( chl_yr$chl < 8.5, 1, 0 )

# Assemble annual TN load data
TNdat <- loads[ which( loads$param=="TN load" & 
                       year(loads$date) >= 2000 # & month(loads$date) %in% 6:9
                       ),
                c('date','value') ]
TNdat$year <- year( TNdat$date )
TN_yr <- TNdat |> group_by(year) |> dplyr::summarise(TN_abs=sum(value)) |> 
  as.data.frame()

# Assemble annual hydrologic load data
hydro_yr <- totanndat[ which( totanndat$bay_segment=="Old Tampa Bay" &
                              totanndat$year >= 2000 ),
                       c('year','hy_load') ] |> as.data.frame()
colnames(hydro_yr) <- c('year','hydro')

# Join data
dat <- inner_join( chl_yr, TN_yr, by = 'year' )
dat <- inner_join( dat, hydro_yr, by = 'year' )

# Compute normalized TN loads
dat$TN_norm <- dat$TN_abs * 449.44 / dat$hydro

# Compute cumulative chl attainment for absolute loads
dat <- dat[ order( dat$TN_abs, decreasing = TRUE ), ]
dat$attain_abs_sum <- c( dat$chl_attain[1], rep(NA,(nrow(dat)-1)) )
for( i in 2:nrow(dat) ){
  dat$attain_abs_sum[i] <- dat$attain_abs_sum[i-1] + dat$chl_attain[i]
}
dat$attain_abs_sum_pct <- dat$attain_abs_sum / nrow(dat)

# Plot attainment ~ TN load (absolute)
png( "../figs/attain_curve_TN_chl.png", width = 7, height = 6, units = 'in', res = 600 )
par(mar=c(5,5.5,3,1))
plot( attain_abs_sum_pct ~ TN_abs, data = dat, xlim = c(300,900),
      type = 'l', lwd = 3, col = rgb(0.1,0.2,0.9,0.8),
      yaxt= 'n', ylab = "", xlab = "TN load (tons/year)",
      main = "Attainment curve (chl < 8.5 ug/L)",
      cex.axis = 1.3, cex.main = 1.3, cex.lab = 1.3 )
axis( 2, at = axTicks(2), labels = paste0(axTicks(2)*100,"%"),
      las = 1, cex.axis = 1.3 )
mtext( "Years attained (%)", side = 2, line = 4, cex = 1.3 )
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
abline( v = seq(0,1000,100), col = rgb(0,0,0,0.1) )
abline( v = seq(0,1000,25), col = rgb(0,0,0,0.1) )
abline( v = 486, col = rgb(0,0,0,0.7), lty = 2, lwd = 2 )
legend( 'topright', bty = 'n', legend = c('absolute','normalized'),
        col = c(rgb(0.1,0.2,0.9,0.8),rgb(1,0.2,0.4,0.8)),
        lwd = 3, cex = 1.3 )

# Compute cumulative chl attainment for normalized loads
dat <- dat[ order( dat$TN_norm, decreasing = TRUE ), ]
dat$attain_norm_sum <- c( dat$chl_attain[1], rep(NA,(nrow(dat)-1)) )
for( i in 2:nrow(dat) ){
  dat$attain_norm_sum[i] <- dat$attain_norm_sum[i-1] + dat$chl_attain[i]
}
dat$attain_norm_sum_pct <- dat$attain_norm_sum / nrow(dat)

# Plot attainment ~ TN load (normalized)
lines( attain_norm_sum_pct ~ TN_norm, data = dat,
       lwd = 3, col = rgb(1,0.2,0.4,0.8) )

dev.off()
