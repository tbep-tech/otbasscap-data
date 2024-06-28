rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
  load( "../data-raw/usgs.RData" )

# Process Lake Tarpon discharge data (USGS)
  usgs$date <- usgs$dateTime |> as.Date()
  # Check continuity over time
  usgs$date |> diff() |> unique()  # there are a few gaps (averaging will take care of it)
  # Check qualifier codes
  usgs$X_00060_00003_cd |> table( useNA = 'always' ) # A & e codes are OK
  usgs$date[ which( usgs$X_00060_00003_cd=="P") ] |> range()  # provisional data in late 2023
  # Calculate monthly average of daily cfs values
  usgs$month <- floor_date( usgs$date, 'month' )
  tarpon <- usgs |> group_by(month) |>
    summarise(discharge_cfs = mean(X_00060_00003)) |> as.data.frame()
  # Convert monthly average cfs to monthly average ft3/day
  tarpon$discharge <- tarpon$discharge_cfs *86400
  # Covert ft3/day to ft3/mo
  tarpon$discharge <- tarpon$discharge * days_in_month(tarpon$month)
  # log transform
  tarpon$logdischarge <- log10( tarpon$discharge+1 )

# Plot monthly discharge data
  par(mfrow=c(2,1))
  plot( discharge ~ month, data = tarpon, type = 'l',
        main = "Lake Tarpon discharge", ylab = "Discharge (ft3/mo)" )
  abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  plot( logdischarge ~ month, data = tarpon, type = 'l',
        main = "Lake Tarpon discharge", ylab = "Discharge (log10 ft3/mo)" )
  abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  
  
# Export monthly discharge data
 save( tarpon, file = "../data-clean/Q_tarpon.RData" )  
 