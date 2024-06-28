rm(list=ls(all=TRUE))

library(dplyr)
library(lubridate)

# load raw data
epcphyto <- readxl::read_excel("../data-raw/PlanktonDataList_ThroughCurrentReportMonth.xlsx")

# subset columns
phyto <- epcphyto |> select( StationNumber, Latitude, Longitude, 
                             SampleTime, LabIDNumber,
                             CLASS, ORDER, FAMILY, NAME,
                             AdjCount, AdjUnits,
                             Qualifier, ID_COMMENTS, qcd
                             ) |> as.data.frame()
colnames( phyto ) <- c( "site","Latitude","Longitude","sampletime","LabID",
                        "class","order","family","name",
                        "cellcount","unit","qualifier","comments","qcd"
                        )

# Correct units
phyto$unit |> unique()
phyto$unit <- "cells/L"

# Add date and day-of-year columns
phyto$date <- phyto$sampletime |> as.Date(format="%Y-%m-%d")
phyto$doy <- phyto$date |> yday()


# Subset data from OTB stations
phyto <- phyto[ which( phyto$site %in% c(36,  38,  40,  41,  42,  46,  47,
                                         50,  51,  60,  61,  62,  63,  64,
                                         65,  66,  67,  68) ), ]

# Label observations by OTB subsegment
phyto$subsegment <- NA
phyto$subsegment[ which( phyto$site %in% c(46,64) ) ] <- "NW"
phyto$subsegment[ which( phyto$site %in% c(47,60,62) ) ] <- "NE"
phyto$subsegment[ which( phyto$site %in% c(42,65,66) ) ] <- "CW"
phyto$subsegment[ which( phyto$site %in% c(40,41,61,63) ) ] <- "CE"
phyto$subsegment[ which( phyto$site %in% c(38,67,68) ) ] <- "SW"
phyto$subsegment[ which( phyto$site %in% c(36,50,51) ) ] <- "SE"

# Remove "TNTC" records (too numerous to count)
phyto$qualifier |> unique()
phyto$cellcount[ which(phyto$qualifier=="TNTC") ] |> table( useNA = 'always' )  # all report zero cells
phyto$name[ which(phyto$qualifier=="TNTC") ] |> table( useNA = 'always' )  # species 
rm.tntc.idx <- which(phyto$qualifier=="TNTC")  
phyto <- phyto[ -rm.tntc.idx, ]  # remove TNTC records

# export
save( phyto, file = "../data-clean/epcphyto.RData" )
