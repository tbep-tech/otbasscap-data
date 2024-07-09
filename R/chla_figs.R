rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

load('../data-clean/epcwq_clean.RData')
epcdat <- epcwq3[which(epcwq3$param == "Chla"),]
epcdat$yr <- year( epcdat$date )
epcdat$month <- floor_date( epcdat$date, 'month' )

subsegs <- c("NW","NE","CW","CE","SW","SE")

# OTB chlorophyll over time -----------------------------------------------

png("../figs/chl_otb.png", res = 600,
    units = 'in', width = 10, height = 8 )
par(mfrow=c(2,1), mar=c(3,4,1,1) )
 # sample data over time
  plot( value ~ date, data = epcdat,
        type = 'l', las = 1,
        col = rgb(0,0.4,0,0.3), lwd = 1,
        main = "",
        ylab = "Chlorophyll a (ug/L)", xlab = '' )
  mtext( "(a) Chlorophyll-a samples across OTB (EPCHC data)",
         side = 3, line = 0, adj = 0, font = 1 )
  abline( v = seq.Date( floor_date(min(epcdat$date),'year'),
                        ceiling_date(max(epcdat$date),'year'), 'year' ),
          col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  points( value ~ date, data = epcdat, cex = 0.4,
          pch = 16, col = rgb(0,0.6,0.4,0.8) )
  points( x = epcdat$date,
          y = rep(-8,nrow(epcdat)),
          pch = "|", cex = 0.4, col = rgb(0.3,0.3,0.3,0.5) )
  # annual mean chlorophyll
  epcmeans <- epcdat |> group_by(yr) |> summarise(value=mean(value)) |>
                as.data.frame()
  plot( value ~ yr, data = epcmeans,
        type = 'l', las = 1,
        col = rgb(0,0.6,0.4,0.8), lwd = 4,
        main = "", ylim = c(5,15), xlim = c(2000,2024),
        ylab = "Chlorophyll a (ug/L)", xlab = '' )
  mtext( "(b) Monthly mean chlorophyll-a concentrations across OTB (EPCHC data)",
         side = 3, line = 0, adj = 0, font = 1 )
  abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  abline( h = c(8.5,9.3), lty = c(3,1), col = rgb(0.3,0.2,0.2,0.7), lwd = 2 )
  points( value ~ yr, data = epcmeans, cex = 1.2, lwd = 2,
          pch = 21, bg = rgb(1,1,1,1), col =  rgb(0,0.6,0.4,1) )
  legend( 'topright', bty = 'n',
          legend = c("Regulatory threshold (9.3 ug/L)",
                     "Management target (8.5 ug/L)"),
          lty = c(1,3), col = rgb(0.3,0.2,0.2,0.7), lwd = 2)
dev.off()


# Sub-segment chlorophyll over time ---------------------------------------

# Plot data and annual averages over time at each OTB subsegment
png("../figs/chl_otbsub_annual.png", res = 600,
    units = 'in', width = 9, height = 6 )
par(mfrow=c(3,2), mar=c(2,4,2,1))
  # Loop through OTB sub-segments
  for( i in 1:length(subsegs) ){
    # Assemble sub-segment data
    this.dat <- epcdat[ which( epcdat$subseg==subsegs[i] ),]
    means <- this.dat |> group_by(yr) |> summarise(value=mean(value)) |> 
               as.data.frame()
    sd <- this.dat |> group_by(yr) |> summarise(sd=sd(value)) |> 
             as.data.frame()
    n <- this.dat |> group_by(yr) |> summarise(n=length(value)) |> 
           as.data.frame()
    means$upr <- means$value + sd$sd/sqrt(n$n)
    means$lwr <- means$value - sd$sd/sqrt(n$n)
    # Generate sub-segment plot
    plot( value ~ yr, data = means,
          type = 'l', las = 1,
          col = rgb(0,0.6,0.4,0.8), lwd = 2,
          main = "", ylim = c(4,25), xlim = c(2000,2024),
          ylab = "Chlorophyll a (ug/L)", xlab = '' )
    mtext( paste0(subsegs[i]," sub-segment"),
           side = 3, line = 0, adj = 0, font = 1 )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
    segments( x0 = means$yr, y0 = means$upr, y1 = means$lwr,
              col = rgb(0,0.6,0.4,1), lwd = 2 )
    abline( h = c(8.5,9.3), lty = c(3,1), col = rgb(0.3,0.2,0.2,0.7), lwd = 2 )
    points( value ~ yr, data = means, cex = 1.2, lwd = 2,
            pch = 21, bg = rgb(1,1,1,1), col =  rgb(0,0.6,0.4,1) )
    if(i==1){
      legend( 'topright', bty = 'n',
              legend = c("Regulatory threshold (9.3 ug/L)",
                         "Management target (8.5 ug/L)"),
              lty = c(1,3), col = rgb(0.3,0.2,0.2,0.7), lwd = 2)
    }
  }
dev.off()


# Sub-segment chlorophyll by month ------------------------------------------------


# Plot data by month-of-year with monthly averages at each OTB subsegment
png("../figs/chl_otbsub_month.png", res = 600,
    units = 'in', width = 9, height = 6 )
par(mfrow=c(3,2), mar=c(2,4,2,1))
  # Loop through OTB sub-segments
  for( i in 1:length(subsegs) ){
    # Assemble sub-segment data
    this.dat <- epcdat[ which( epcdat$subseg==subsegs[i] ),]
    this.dat$month <- this.dat$date |> month()
    means <- this.dat |> group_by(month) |> summarise(value=mean(value)) |> 
               as.data.frame()
    sd <- this.dat |> group_by(month) |> summarise(sd=sd(value)) |> 
            as.data.frame()
    n <- this.dat |> group_by(month) |> summarise(n=length(value)) |> 
            as.data.frame()
    means$upr <- means$value + sd$sd/sqrt(n$n)
    means$lwr <- means$value - sd$sd/sqrt(n$n)
    # Generate sub-segment plot
    plot( value ~ month, data = means,
          type = 'l', las = 1,
          col = rgb(0,0.6,0.4,0.8), lwd = 2,
          main = "", ylim = c(0,35), xaxt = 'n',
          ylab = "Chlorophyll a (ug/L)", xlab = '' )
    mtext( paste0(subsegs[i]," sub-segment  (",
                  min(this.dat$yr),"\u2012",max(this.dat$yr),")"),
           side = 3, line = 0, adj = 0, font = 1 )
    axis( 1, at = means$month, labels = month(means$month,label=TRUE) )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
    segments( x0 = means$month, y0 = means$upr, y1 = means$lwr,
              col = rgb(0,0.6,0.4,1), lwd = 2 )
    abline( h = c(8.5,9.3), lty = c(3,1), col = rgb(0.3,0.2,0.2,0.7), lwd = 2 )
    points( value ~ month, data = means, cex = 1.2, lwd = 2,
            pch = 21, bg = rgb(1,1,1,1), col =  rgb(0,0.6,0.4,1) )
    if(i==1){
      legend( 'topright', bty = 'n',
              legend = c("Regulatory threshold (9.3 ug/L)",
                         "Management target (8.5 ug/L)"),
              lty = c(1,3), col = rgb(0.3,0.2,0.2,0.7), lwd = 2)
    }
  } # // end i loop
dev.off()



