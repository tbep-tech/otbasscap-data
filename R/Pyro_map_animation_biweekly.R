rm(list=ls(all=TRUE)) 

library(lubridate)
library(rgdal)
library(sf)
library(animation)


# Load Pyro data from file
load( "../data/Pyro.Rdata")
TB.map <- readOGR( dsn = "../data/Tampa_Bay_Shoreline",
                   layer = "Tampa_Bay_Shoreline" )

# Subset and process Pyro data
  pyro <- pyro[ which( pyro$yr >= 2012 ), ]
  pyro$logcells <- log10( pyro$pyro )
  pyro$logcells[which(is.na(pyro$logcells))] <- 0
  # Define function to floor dates to 'biweekly' time frame
  biweekly <- function( d ){
    library(lubridate)
    if( mday(d)<15 ){
      mo <- month(d)
      if(mo<10){
        mo <- paste0("0",mo)
      }
      out <- paste0( year(d),"-",mo,"-01")
    } else {
      mo <- month(d)
      if(mo<10){
        mo <- paste0("0",mo)
      }
      out <- paste0( year(d),"-",mo,"-15")
    }
    return(out)
  }  # // end biweekly()
  # Apply biweekly function
  pyro$week <- apply( data.frame(pyro$date), 1, biweekly ) |> as.Date()


# Create animated map file (mp4) by looping through weeks
# to generate animation frames
weeks <- seq.Date( floor_date(min(pyro$date),'month'),
                   ceiling_date(max(pyro$date),'month'),
                   'day'
                  ) |> as.data.frame() |> apply(1,biweekly) |> unique()

# Download ffmpeg at https://ffmpeg.org/download.html
ani.options(ffmpeg = "C:/Users/miles/Documents/GitHub/otbasscap-data/R/ffmpeg.exe")

t.start <- Sys.time()
t.start

saveVideo(
  expr = {
    ani.record(reset = TRUE) # discard previous plots
    ani.options(interval = 1/2)
  
    for( i in 1:length(weeks) ){
      # subset data for this week
      this <- pyro[ which( pyro$week==weeks[i]), ]
      # draw map
      plot( TB.map, border = rgb(0,0,0,0.3), # col = rgb(0.2,0.2,0.4,0.3),
            xlim = c(-82.75,-82.50), ylim = c(27.80,28.05),
            main = paste0("Pyrodinium sampling at Old Tampa Bay (biweekly bins)\n",
                          weeks[i],
                          # month(weeks[i],label=TRUE,abbr=TRUE),
                          # ' ',year(weeks[i]),
                          '\nN = ',nrow(this) ) )
      # plot sampling points
      if( length(this)>0 ){
        points( Latitude ~ Longitude,  # zero samples
                data = this[which(this$logcells==0),],
                pch = 1, col = rgb(1,0.1,0.1,0.6), lwd = 1,
                cex = 1 )
        points( Latitude ~ Longitude,  # non-zero samples
                data = this[which(this$logcells>0),],
                pch = 16, col = rgb(1,0.1,0.1,0.6),
                cex = this$logcells[which(this$logcells>0)] )
      }
      # legend
      legend( 'bottom', horiz = TRUE,
              title = "Abundance (cells/L)",
              title.font = 2,
              legend = c(0,expression(10^2~"     "),
                           expression(10^3), expression(10^4),
                           expression(10^5), expression(10^6)),
              pch = c(1,rep(16,5)), pt.cex = c(1,2:6),
              col = rgb(1,0.1,0.1,0.6)
              )
    }  # // end i loop
    
  },
  
  video.name = "../figs/Pyro_map_animation_biweekly.mp4",
  ani.width = 3600,
  ani.height = 3600,
  ani.res = 400

)  # // end saveVideo

t.end <- Sys.time()
t.end - t.start
