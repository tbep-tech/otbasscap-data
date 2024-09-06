rm(list=ls(all=TRUE)) 

library(lubridate)

# Read FWC Pyrodinium data from file
pyro <- read.csv("../data-raw/fwcpyro20112023.csv")

# Coerce date column to Date
pyro$date <- pyro$date |> as.Date(format="%Y-%m-%d")

# Label observations by OTB subsegment
pyro$subsegment <- NA
pyro$subsegment[ which( pyro$Latitude > 27.97 & pyro$Longitude < -82.66 ) ] <- 'NW'
pyro$subsegment[ which( pyro$Latitude > 27.969 & pyro$Longitude > -82.64 ) ] <- 'NE'
pyro$subsegment[ which( pyro$Latitude < 27.95 & pyro$Latitude > 27.92 &
                        pyro$Longitude < -82.62 ) ] <- 'CW'
pyro$subsegment[ which( (pyro$Latitude > 27.95  & pyro$Latitude < 27.97 &
                         pyro$Longitude > -82.63 & pyro$Longitude < -82.57) |
                        (pyro$Latitude > 27.94  & pyro$Latitude < 27.95 &
                         pyro$Longitude > -82.58 & pyro$Longitude < -82.56) ) ] <- 'CE'
pyro$subsegment[ which( pyro$Latitude > 27.84 & pyro$Latitude < 27.91 ) ] <- 'SW'
pyro$subsegment[ which( pyro$Longitude > -82.55 ) ] <- 'SE'

# Remove records outside of OTB (is.na(pyro$subsegment)==TRUE)
pyro <- pyro[ -which( is.na(pyro$subsegment)), ]

# Label records as routine (TRUE) or event-based (FALSE) samples
  # concatenate lat-lon coords into strings
  coordstr <- paste( pyro$Latitude, pyro$Longitude, sep = ", " )
  # tabulate coord strings; top 9 are the routine sites (confirmed independently)
  routine_sites <- ( table(coordstr) |> sort(decreasing=TRUE) |> names() )[1:9]
  # create new column to label routine samples (TRUE)
  pyro$routine <- FALSE
  pyro$routine[ which( coordstr %in% routine_sites ) ] <- TRUE
  # plot to confirm (routine locations in red)
  plot( Latitude ~ Longitude, data = pyro, cex = 2 )
  points( Latitude ~ Longitude, data = pyro[which(pyro$routine==TRUE),],
          cex = 2, col = rgb(1,0,0,1), pch=16 )

# export as RData and csv
save( pyro, file = "../data/pyro.RData" )
write.csv( pyro, file = "../data/pyro.csv", row.names = FALSE )