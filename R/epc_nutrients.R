rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# Load & process chl data
load('../data-clean/epcwq_clean.RData')
epcdat <- epcwq3[which(epcwq3$param %in% c("TN","NOx","TP","PO4") ),]
epcdat$yr <- year( epcdat$date )
epcdat$month <- floor_date( epcdat$date, 'month' )

# Define geometric mean function
geomean <- function(x){
  return( 10^mean(log10(x)) )
}

# Specify OTB sub-segments
subsegs <- c("NW","NE","CW","CE","SW","SE")

# Sub-segment nitrogen conc over time ---------------------------------------

# Plot data and annual averages over time at each OTB subsegment
png("../figs/epc_nitrogen_otbsub_annual.png", res = 600,
    units = 'in', width = 9, height = 6 )
par(mfrow=c(3,2), mar=c(2,4,2,1))
  # Loop through OTB sub-segments
  for( i in 1:length(subsegs) ){
    # Assemble sub-segment TN data
    this.TN.dat <- epcdat[ which( epcdat$subseg==subsegs[i] & epcdat$param=="TN" ),]
    TN.means <- this.TN.dat |> group_by(yr) |> dplyr::summarise(value=geomean(value)) |> 
               as.data.frame()
    TN.sd <- this.TN.dat |> group_by(yr) |> dplyr::summarise(sd=sd(value)) |> 
             as.data.frame()
    TN.n <- this.TN.dat |> group_by(yr) |> dplyr::summarise(n=length(value)) |> 
           as.data.frame()
    TN.means$upr <- TN.means$value + TN.sd$sd/sqrt(TN.n$n)
    TN.means$lwr <- TN.means$value - TN.sd$sd/sqrt(TN.n$n)
    # Assemble sub-segment NOx data
    this.NOx.dat <- epcdat[ which( epcdat$subseg==subsegs[i] & epcdat$param=="NOx" ),]
    NOx.means <- this.NOx.dat |> group_by(yr) |> dplyr::summarise(value=geomean(value)) |> 
      as.data.frame()
    NOx.sd <- this.NOx.dat |> group_by(yr) |> dplyr::summarise(sd=sd(value)) |> 
      as.data.frame()
    NOx.n <- this.NOx.dat |> group_by(yr) |> dplyr::summarise(n=length(value)) |> 
      as.data.frame()
    NOx.means$upr <- NOx.means$value + NOx.sd$sd/sqrt(NOx.n$n)
    NOx.means$lwr <- NOx.means$value - NOx.sd$sd/sqrt(NOx.n$n)
    # Generate sub-segment plot
    # TN (black)
    plot( value ~ yr, data = TN.means,
          type = 'l', las = 1,
          col = rgb(0,0,0,0.8), lwd = 2,
          main = "", ylim = c(0,1.25), xlim = c(2000,2024),
          ylab = "Concentration (mg/L)", xlab = '' )
    mtext( paste0(subsegs[i]," sub-segment"),
           side = 3, line = 0, adj = 0, font = 1 )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
    segments( x0 = TN.means$yr, y0 = TN.means$upr, y1 = TN.means$lwr,
              col = rgb(0,0,0,1), lwd = 2 )
    points( value ~ yr, data = TN.means, cex = 1.2, lwd = 2,
            pch = 21, bg = rgb(1,1,1,1), col =  rgb(0,0,0,1) )
    # NOx
    lines( value ~ yr, data = NOx.means, lwd = 2, col = rgb(1,0.2,0.1,0.8) )
    segments( x0 = NOx.means$yr, y0 = NOx.means$upr, y1 = NOx.means$lwr,
              col = rgb(1,0.2,0.1,0.8), lwd = 2 )
    points( value ~ yr, data = NOx.means, cex = 1.2, lwd = 2,
            pch = 21, bg = rgb(1,1,1,1), col =  rgb(1,0.2,0.1,1) )
    if(i==1){
      legend( 'topright', bty = 'n',
              legend = c("TN","NOx"),
              col = c( 1, rgb(1,0.2,0.1,1) ),
              lwd = 2 )
    }
  }
dev.off()


# Sub-segment phosphorus over time ---------------------------------------

# Plot data and annual averages over time at each OTB subsegment
png("../figs/epc_phosphorus_otbsub_annual.png", res = 600,
    units = 'in', width = 9, height = 6 )
par(mfrow=c(3,2), mar=c(2,4,2,1))
# Loop through OTB sub-segments
for( i in 1:length(subsegs) ){
  # Assemble sub-segment TP data
  this.TP.dat <- epcdat[ which( epcdat$subseg==subsegs[i] & epcdat$param=="TP" ),]
  TP.means <- this.TP.dat |> group_by(yr) |> dplyr::summarise(value=geomean(value)) |> 
    as.data.frame()
  TP.sd <- this.TP.dat |> group_by(yr) |> dplyr::summarise(sd=sd(value)) |> 
    as.data.frame()
  TP.n <- this.TP.dat |> group_by(yr) |> dplyr::summarise(n=length(value)) |> 
    as.data.frame()
  TP.means$upr <- TP.means$value + TP.sd$sd/sqrt(TP.n$n)
  TP.means$lwr <- TP.means$value - TP.sd$sd/sqrt(TP.n$n)
  # Assemble sub-segment PO4 data
  this.PO4.dat <- epcdat[ which( epcdat$subseg==subsegs[i] & epcdat$param=="PO4" ),]
  PO4.means <- this.PO4.dat |> group_by(yr) |> dplyr::summarise(value=geomean(value)) |> 
    as.data.frame()
  PO4.sd <- this.PO4.dat |> group_by(yr) |> dplyr::summarise(sd=sd(value)) |> 
    as.data.frame()
  PO4.n <- this.PO4.dat |> group_by(yr) |> dplyr::summarise(n=length(value)) |> 
    as.data.frame()
  PO4.means$upr <- PO4.means$value + PO4.sd$sd/sqrt(PO4.n$n)
  PO4.means$lwr <- PO4.means$value - PO4.sd$sd/sqrt(PO4.n$n)
  # Generate sub-segment plot
  # TP (black)
  plot( value ~ yr, data = TP.means,
        type = 'l', las = 1,
        col = rgb(0,0,0,0.8), lwd = 2,
        main = "", ylim = c(0,0.25), xlim = c(2000,2024),
        ylab = "Concentration (mg/L)", xlab = '' )
  mtext( paste0(subsegs[i]," sub-segment"),
         side = 3, line = 0, adj = 0, font = 1 )
  abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  segments( x0 = TP.means$yr, y0 = TP.means$upr, y1 = TP.means$lwr,
            col = rgb(0,0,0,1), lwd = 2 )
  points( value ~ yr, data = TP.means, cex = 1.2, lwd = 2,
          pch = 21, bg = rgb(1,1,1,1), col =  rgb(0,0,0,1) )
  # PO4
  lines( value ~ yr, data = PO4.means, lwd = 2, col = rgb(1,0.2,0.1,0.8) )
  segments( x0 = PO4.means$yr, y0 = PO4.means$upr, y1 = PO4.means$lwr,
            col = rgb(1,0.2,0.1,0.8), lwd = 2 )
  points( value ~ yr, data = PO4.means, cex = 1.2, lwd = 2,
          pch = 21, bg = rgb(1,1,1,1), col =  rgb(1,0.2,0.1,1) )
  if(i==1){
    legend( 'topright', bty = 'n',
            legend = c("TP","PO4"),
            col = c( 1, rgb(1,0.2,0.1,1) ),
            lwd = 2 )
  }
}
dev.off()