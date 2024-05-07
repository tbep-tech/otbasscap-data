rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# load data
load( here::here("data/epcwq.RData"))
epcwq1 <- epcwq |> as.data.frame()

# format dates
epcwq1$date <- epcwq1$SampleTime |> as.Date("%Y-%m-%d")

# subset columns for selected wq parameters
# --> data and 'Q' columns must be listed in consecutive pairs!!
epcwq1 <- select( epcwq1, date, StationNumber, SampleDepth,
                  Total_Nitrogen, Total_NitrogenQ,
                  Chlorophyll_a, Chlorophyll_aQ,
                  `DO-T`, `DOQ-T`,
                  `DO-M`, `DOQ-M`,
                  `DO-B`, `DOQ-B`,
                  `Sal-T`, `SalQ-T`,
                  `Sal-M`, `SalQ-M`,
                  `Sal-B`, `SalQ-B`,
                  SecchiDepth, Secchi_Q,
                  Total_Suspended_Solids, Total_Suspended_SolidsQ,
                  Turbidity, TurbidityQ,
                  `TempWater-T`, `TempWaterQ-T`,
                  `TempWater-M`, `TempWaterQ-M`,
                  `TempWater-B`, `TempWaterQ-B` )

# subset OTB stations
epcwq1 <- epcwq1[ which( epcwq1$StationNumber %in% c(36,  38,  40,  41,  42,  46,  47,
                                                     50,  51,  60,  61,  62,  63,  64,
                                                     65,  66,  67,  68) ), ]

# subset by date
epcwq1 <- epcwq1[ which( year(epcwq1$date) >= 2000 ), ]

# loop through data and QC columns to pivot from wide to long format
datacols <- seq( 4, 30, 2 )
epcwq2 <- data.frame( matrix(ncol=7,nrow=0) )
for( i in datacols ){
  temp <- epcwq1[ ,1:3 ]
  temp$Analyte <- colnames(epcwq1)[i]
  temp$value <- epcwq1[,i] |> as.numeric()
  temp$unit <- NA
  temp$QA <- epcwq1[,i+1]
  temp <- temp[ complete.cases( temp$value ), ]
  epcwq2 <- rbind( epcwq2, temp )
}

# standardize parameter names
  param.names.old <- epcwq2$Analyte |> unique() |> sort()
  param.names.new <- c( "Chla", "DO_bot", "DO_mid", "DO_top",
                        "Sal_bot", "Sal_mid", "Sal_top","Secchi",
                        "Temp_bot","Temp_mid","Temp_top",
                        "TN", "TSS", "Turbidity" )
  epcwq2$param <- mapvalues( epcwq2$Analyte, param.names.old, param.names.new )
  
  # assign units
  units <- c( "ug/l", "mg/l", "mg/l", "mg/l",
              "PSU", "PSU", "PSU", "m",
              "ppt", "ppt", "ppt",
              "mg/l", "mg/l", "NTU" )  # listed in same order as param.names.new
  epcwq2$unit <- mapvalues( epcwq2$param, param.names.new, units )

# qualifier codes
  # replace "NULL" with NA
  epcwq2$QA[ which(epcwq2$QA=="NULL") ] <- NA
  # specify codes
  fatal.codes <- c("*","?","A","B","G","H","J","K",
                   "L","N","O","Q","T","V","Y","Z")
  # define function to locate fatal records
  find.fatal <- function( QUALIFIER, FATAL=fatal.codes ){  # QUALIFER arg is a string
    input.code <- QUALIFIER |> strsplit(split='') |>
      unlist() |> toupper()  # parse string into single characters
    fatal <- input.code %in% FATAL |> any()  # check if any characters match FATAL
    return( fatal )  # return TRUE or FALSE
  }  # // end find.fatal()
  # apply function to locate fatal records
  rm.fatal.idx <- apply( matrix(epcwq2$QA,ncol=1), 1, find.fatal ) |> which()
  epcwq2 <- epcwq2[ -rm.fatal.idx, ]

# acknowledge Secchi records visible on bottom
  epcwq2[ which( epcwq2$param=="Secchi" & epcwq2$QA==">"), ] |> nrow() # number of VOB Secchi records
  epcwq2[ which( epcwq2$param=="Secchi"), ] |> nrow()  # total number of Secchi records

# coerce depth to numeric
  epcwq2$SampleDepth <- epcwq2$SampleDepth |> as.numeric()

# check for duplicates (confirm FALSE)
  epcwq2 |> duplicated() |> any()
  
# check for non-positive values (confirm FALSE)
  any( epcwq2$value == 0 )

# subset and rename columns
  epcwq2 <- epcwq2[, c("date","param","value","unit","StationNumber") ]
  colnames( epcwq2 ) <- c("date","param", "value", "unit", "site" )
  
# aggregate data to monthly timeframe
  epcwq3 <- epcwq2
  epcwq3$date <- epcwq3$date |> floor_date('month')
  epcwq3 <- epcwq2 |> dplyr::summarise( value = mean(value),
                                        .by = c(date,param,site,unit) )

# export clean dataset
  save( epcwq3, file = "data-clean/epcwq_clean.RData" )

  