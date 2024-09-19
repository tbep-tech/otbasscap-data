# Alternative reference periods
rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# Load hydro & TN load data
load("../data-raw/totanndat.RData")
loads <- totanndat[ which( totanndat$bay_segment=="Old Tampa Bay"),
                              c('year','hy_load','tn_load') ] |> as.data.frame()

# Calculate alternative delivery ratios
  # Specify reference periods (start and end years)
  ref_periods <- data.frame( year_str = c( 1992, 2000, 2004, 2004, 2007, 2015 ),
                             year_end = c( 1994, 2021, 2010, 2006, 2009, 2021 ),
                             comment = c("Current paradigm",
                                         "Broad contemparary period",
                                         "Precedes historic seagrass growth/recovery",
                                         "3-year subset of 2004-2010 seagrass growth period",
                                         "3-year subset of 2004-2010 seagrass growth period",
                                         "Precedes/includes seagrass decline")
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
