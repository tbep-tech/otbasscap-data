
# This script queries USGS rainfall data at the site closest to OTB
# and estimates total annual rainfall. The purpose was to estimate
# hydrologic conditions in 2022 and 2023 compared to other years,
# since hydrologic load estimates from TBEP are currently only avaiable 
# through 2021.

rm(list=ls(all=TRUE)) 

library( dataRetrieval )
library( lubridate )
library( dplyr )

# Query USGS data
dat.usgs <- readNWISdata( site = "280842082392000",
                          parameterCd = c( "00045" ),
                          startDate = "2000-01-01",
                          endDate = "2023-12-31"
                          )

# Process data
dat <- data.frame( date = as.Date(dat.usgs$dateTime),
                   rain = dat.usgs$X_00045_00006,
                   cd = dat.usgs$X_00045_00006_cd
                   )
range( dat$date )

dat$month <- floor_date( dat$date )
dat$year <- year( dat$date )
table( dat$year )  # many years do not have complete daily data

# Compute each year's daily average and multiply by 365 days to estimate annual totals
rain <- dat |> group_by(year) |> dplyr::summarise(dailymean=mean(rain)) |> as.data.frame()
rain$anntot <- rain$dailymean * 365

# Generate Plot
plot( anntot ~ year, data = rain, type = 'l', lwd = 3,
      main = "Annual rainfall at Brooker Creek near Tarpon Springs (USGS 280842082392000)",
      xlab = "", ylab = "Total rainfall (inches)", las = 1 )
points( anntot ~ year, data = rain, pch = 21,
        bg = rgb(1,1,1,1,), lwd = 2 )
abline( v = 2000:2024, col = rgb(0,0,0,0.1) )
abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
