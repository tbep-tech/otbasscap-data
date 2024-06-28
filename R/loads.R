# This script aggregates raw monthly load estimates for OTB (by source)
# to produce total monthly TN, TP, TSS, and BOD loads across all sources
# 2024, Miles Medina, ECCO Scientific
rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# Load data
loads.raw <- read.csv( "../data-raw/TB_loads_monthly.csv" )
loads.raw <- loads.raw[ which(loads.raw$bay_segment=="Old Tampa Bay"), ]

# Rename column
colnames(loads.raw)[5:8] <- sub("_load","",colnames(loads.raw)[5:8])
colnames(loads.raw)[5:8] <- colnames(loads.raw)[5:8] |> 
                              toupper() |> unlist() |> paste("load")

# Loop through load constituents to assemble total monthly loads
loads <- data.frame( date = Date(),
                     value = as.numeric(),
                     param = as.character(),
                     site = as.character(),
                     unit = as.character()
                     ) # initiate loads dataframe
vars <- colnames(loads.raw)[5:8]
for( i in 1:length(vars) ){
  
  # subset data for i'th contituent
  this <- loads.raw[ ,c('year','month',vars[i]) ]
  colnames(this)[3] <- "load"
  # assemble Date
  this$date <- as.Date( paste0(this$year,"-",this$month,"-",01), format="%Y-%m-%d" )
  # aggregate loads across source types
  this <- this |> group_by(date) |> summarise(value=sum(load)) |> as.data.frame()
  # create labels
  this$param <- vars[i]
  this$site <- "OTB watershed"
  this$unit <- "tons"
  # append this to loads dataframe
  loads <- rbind( loads, this )
  
}  # // end vars loop

# Plot monthly loads
par(mfrow=c(4,1), mar=c(3,4,2,1))
for( i in 1:length(vars) ){
  plot( value ~ date, data = loads[which(loads$param==vars[i]),],
        las = 1, type = 'l', lwd = 2, col = rgb(0,0,0,0.7),
        main = vars[i], xlab = "", ylab = "load (tons)"
        )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  abline( v = seq.Date( min(loads$date),
                        max(loads$date)+60, 'year'),
          col = rgb(0,0,0,0.1) )
  
}  # end plotting loop

# Export load data
save( loads, file = "../data/loads.RData" )