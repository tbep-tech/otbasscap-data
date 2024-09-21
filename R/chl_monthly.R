rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

png( "../figs/chl_monthly.png", height = 10, width = 10, units = 'in', res = 600 )
par(mfrow=c(2,1),mar=c(4,4,1,1))

# Monthly mean chl time series for OTB ------------------------------------
load( "../data-clean/epcwq_clean.RData" )
dat <- epcwq3[ which( epcwq3$param=="Chla" ), ]
chl_mo <- dat[,c('date','value')] |> group_by(date) |> 
  dplyr::summarise(value=mean(value)) |> as.data.frame()
attain <- c(2000,2002,2005,2006,2010,2012,2014,2022,2023)  # chl attainment years
plot( value ~ date, data = chl_mo,
      type = 'l', ylim = c(0,40),
      xlab = "", ylab = "Chlorophyll a (ug/L)",
      col = rgb(0.1,0.7,0.4,0.7), lwd = 2,
      las = 1, cex.axis = 1.3, cex.lab = 1.3
)
for( i in 1:length(attain) ){
  polygon( x = c( as.Date(c(paste0(attain[i],"-01-01"),
                            paste0(attain[i],"-12-31"),
                            paste0(attain[i],"-12-31"),
                            paste0(attain[i],"-01-01"))) ),
           y = c(0,0,50,50), col = rgb(0,0.6,0.8,0.1), border = rgb(0,0,0,0) )
}
mtext( "(a) Monthly mean chlorophyll-a concentration at OTB", side = 3, line = 0,
       font = 1, cex = 1.3, adj = 0 )
# abline( v = seq.Date(as.Date("2000-01-01"), as.Date("2024-01-01"), "6 months"), col=rgb(0,0,0,0.1) )
abline( v = seq.Date(as.Date("2000-01-01"), as.Date("2024-01-01"), "year"), col=rgb(0,0,0,0.2) )
abline( h = seq(0,50,5), col = rgb(0,0,0,0.1))
abline( h=0 )


# Chlorophyll signal (SSA) by month ---------------------------------------
load("../data/ssa_chl_otb.RData")
ssa_chl_otb$dat$mo <- ssa_chl_otb$dat$month |> month()
colors <- hcl.colors( 24, "viridis", alpha = 0.7 )
  # Initialize sub-segment plot
  plot( signal ~ mo, dat = ssa_chl_otb$dat,
        type = 'l', ylim = c(0,40),
        las = 1, cex.axis = 1.3, cex.lab = 1.3,
        xlab = "", ylab = "Chlorophyll a (ug/L)", xaxt = 'n',
        lwd = 2, col = rgb(0,0,0,0) )
  axis( 1, at = 1:12, labels = month(1:12,label=TRUE), cex.axis = 1.3 )
  abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  mtext( "(b) Chlorophyll-a concentration signal at OTB, by month", side = 3, line = 0,
         font = 1, cex = 1.3, adj = 0 )
  # Loop through years to plot subseg signal by month
  yrs <- ssa_chl_otb$dat$month |> year() |> unique() |> sort()
  for( j in 1:length(yrs) ){
    this.j <- ssa_chl_otb$dat[ which( year(ssa_chl_otb$dat$month)==yrs[j] ), ]
    lines( signal ~ mo, data = this.j,
           lwd = 3, col = colors[j] )
  }  # // end j loop
  # Draw legend
    polygon( x = c(1,2,2,1), y = c(40,40,30,30),
             col = rgb(1,1,1,1), border = rgb(0.3,0.3,0.3,1) )
    as.raster(colors,ncol=1) |> rasterImage(1,30,2,40)
    text( x = 2, y = 39.5, label = "2000", pos = 4 )
    text( x = 2, y = 30.5, label = "2023", pos = 4 )
    
dev.off()


# Monthly mean chl time series for OTB sub-segments ----------------------------
load( "../data-clean/epcwq_clean.RData" )
dat <- epcwq3[ which( epcwq3$param=="Chla" ), ]
subsegs <- c("NW","NE","CW","CE","SW","SE")
png( "../figs/chl_monthly_subseg.png", height = 10, width = 10, units = 'in', res = 600 )
par(mfrow=c(3,2),mar=c(4,4,2,2))
for( subseg in subsegs ){
  
  this.dat <- dat[ which(dat$subseg==subseg), ]
  
  this.chl <- this.dat[,c('date','value')] |> group_by(date) |> 
    dplyr::summarise(value=mean(value)) |> as.data.frame()
  
  plot( value ~ date, data = this.chl,
        type = 'l', ylim = c(0,100),
        xlab = "", ylab = "Chlorophyll a (ug/L)",
        col = rgb(0.1,0.7,0.4,0.7), lwd = 2,
        las = 1, cex.axis = 1.3, cex.lab = 1.3
  )
  for( i in 1:length(attain) ){
    polygon( x = c( as.Date(c(paste0(attain[i],"-01-01"),
                              paste0(attain[i],"-12-31"),
                              paste0(attain[i],"-12-31"),
                              paste0(attain[i],"-01-01"))) ),
             y = c(0,0,200,200), col = rgb(0,0.6,0.8,0.1), border = rgb(0,0,0,0) )
  }
  mtext( paste0("Monthly mean chlorophyll-a (",subseg," sub-segment)"),
         side = 3, line = 0,
         font = 1, cex = 1, adj = 0 )
  # abline( v = seq.Date(as.Date("2000-01-01"), as.Date("2024-01-01"), "6 months"), col=rgb(0,0,0,0.1) )
  abline( v = seq.Date(as.Date("2000-01-01"), as.Date("2024-01-01"), "year"), col=rgb(0,0,0,0.2) )
  abline( h = seq(0,200,10), col = rgb(0,0,0,0.1))
  abline( h=0 )
}
dev.off()


# Chlorophyll attainment and Pyro -----------------------------------------

load( "../data/Pyro.Rdata")
pyro <- pyro[ which(pyro$routine==TRUE & year(pyro$date)>=2012), c('date','pyro') ]
pyro$pyro[ which(is.na(pyro$pyro)) ] <- 0  # NA means 0 cells
pyro$month <- floor_date( pyro$date, 'month' )
pyro_med <- pyro |> group_by(month) |> dplyr::summarise(value=median(pyro)) |> as.data.frame()
pyro_max <- pyro |> group_by(month) |> dplyr::summarise(value=max(pyro)) |> as.data.frame()
pyro_90p <- pyro |> group_by(month) |> dplyr::summarise(value=quantile(pyro,0.90)) |> as.data.frame()
pyro_med$value <- log10(pyro_med$value+1)
pyro_max$value <- log10(pyro_max$value+1)
pyro_90p$value <- log10(pyro_90p$value+1)
png( "../figs/chl_monthly_pyro.png", height = 4, width = 8, units = 'in', res = 600 )
par(mfrow=c(1,1),mar=c(4,4,2,2))
plot( value ~ month, data = pyro_max,
      type = 'l', ylim = c(0,7),
      xlim = c(as.Date("2012-01-01"),as.Date("2024-01-01")),
      xlab = "", ylab = "P. bahamense (cells/L)",
      col = rgb(0,0,0,0), lwd = 2, xaxt = 'n', yaxt = 'n',
      las = 1, cex.axis = 1.3, cex.lab = 1.3
)
axis( 1, at = seq.Date(as.Date("2012-01-01"),as.Date("2024-01-01"),'year'),
      labels = 2012:2024, cex.axis = 1.3 )
axis( 2, at = 0:7, las = 1, cex.axis = 1.3,
      labels = c( 0, expression("10"^1), expression("10"^2),
                  expression("10"^3), expression("10"^4),
                  expression("10"^5), expression("10"^6), 
                  expression("10"^7) ) )
for( i in 1:length(attain) ){
  polygon( x = c( as.Date(c(paste0(attain[i],"-01-01"),
                            paste0(attain[i],"-12-31"),
                            paste0(attain[i],"-12-31"),
                            paste0(attain[i],"-01-01"))) ),
           y = c(0,0,50,50), col = rgb(0,0.6,0.8,0.1), border = rgb(0,0,0,0) )
}
# abline( v = seq.Date(as.Date("2000-01-01"), as.Date("2024-01-01"), "6 months"), col=rgb(0,0,0,0.1) )
abline( v = seq.Date(as.Date("2000-01-01"), as.Date("2024-01-01"), "year"), col=rgb(0,0,0,0.2) )
abline( h = 0:10, col = rgb(0,0,0,0.1))
polygon( x = c( pyro_max$month, rep(0,nrow(pyro_max)) ),
         y = c( pyro_max$value, rev( pyro_max$value ) ),
         col = rgb(1,0.4,0.5,1), border = rgb(0,0,0,0) )
polygon( x = c( pyro_med$month, rep(0,nrow(pyro_med)) ),
         y = c( pyro_med$value, rev( pyro_med$value ) ),
         col = rgb(1,0.6,0.7,1), border = rgb(1,1,1,0.4) )

# lines( value ~ month, data = pyro_max, type = 'l', col = rgb(1,0.2,0.1,0.7), lwd = 2 )
# lines( value ~ month, data = pyro_90p, type = 'l', col = rgb(1,0.4,0.2,0.7), lwd = 2 )
# lines( value ~ month, data = pyro_med, type = 'l', col = rgb(0.6,0.3,0.3,0.7), lwd = 2 )
# mtext( "Monthly max P. bahamense abundance w/ chl attainment (8.5 ug/l)", side = 3, line = 0,
#        font = 1, cex = 1.3, adj = 0 )
abline( h=0 )
legend( 'top', horiz = TRUE, bty = 'n',
        legend = c("Maximum cell count   ","Median cell count"),
        fill = c(rgb(1,0.4,0.5,1),rgb(1,0.6,0.7,1)) 
       )
dev.off()

