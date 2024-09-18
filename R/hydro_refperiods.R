# Alternative reference periods
rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# Load hydro & TN load data
load("../data-raw/totanndat.RData")
loads <- totanndat[ which( totanndat$bay_segment=="Old Tampa Bay"),
                              c('year','hy_load','tn_load') ] |> as.data.frame()

# # Load data
# hydro <- read.csv("../data-raw/TB_hydro_monthly.csv")
# load('../data/loads.RData')
# 
# # Assemble annual hydrologic loads
# hydro$date <- as.Date( paste0(hydro$year,"-",formatC(hydro$month,width=2,flag="0"),"-01") )
# hydro <- hydro[ which(hydro$bay_segment=="Old Tampa Bay" ), ]
# hydro <- hydro[, c("date","hy_load_106_m3_mo") ]
# hydro$year <- year( hydro$date )
# hydro <- hydro |> group_by(year) |> dplyr::summarise(hy_load_106_m3_yr=sum(hy_load_106_m3_mo)) |> 
#            as.data.frame()
# plot( hy_load_106_m3_yr ~ year, data = hydro, type = 'b' )
# 
# # Assemble annual TN loads
# loads <- loads[ which( loads$param=="TN load" ), ]
# loads$year <- year( loads$date )
# TN_load <- loads |> group_by(year) |> dplyr::summarise(tn_load_tons_yr=sum(value)) |> 
#   as.data.frame()


# Calculate alternative delivery ratios
  # Specify reference periods (start and end years)
  ref_periods <- data.frame( year_str = c( 1992, 2000, 2004, 2015 ),
                             year_end = c( 1994, 2021, 2010, 2021 ),
                             comment = c("current paradigm",
                                         "broad contemparary period",
                                         "precedes seagrass growth/recovery",
                                         "precedes/includes seagrass decline")
                            )
  # Initialize delivery ratio table
  deliv <- data.frame( year_str = ref_periods$year_str,
                       year_end = ref_periods$year_end,
                       TN_mean = rep(NA,nrow(ref_periods)),
                       hydro_mean = rep(NA,nrow(ref_periods)),
                       delivery_ratio = rep(NA,nrow(ref_periods)),
                       comment = ref_periods$comment
                       )
  # Populate delivery ratio table
  for( i in 1:nrow(deliv) ){
    deliv$TN_mean[i] <- mean( loads$tn_load[ which(loads$year %in% deliv$year_str[i]:deliv$year_end[i] ) ] )
    deliv$hydro_mean[i] <- mean( loads$hy_load[ which(loads$year %in% deliv$year_str[i]:deliv$year_end[i] ) ] )
    deliv$delivery_ratio[i] <- deliv$TN_mean[i] / deliv$hydro_mean[i]
  }
  
# Export
  deliv
  write.csv( deliv, file = "../data/hydro_refperiods.csv", row.names = FALSE )
