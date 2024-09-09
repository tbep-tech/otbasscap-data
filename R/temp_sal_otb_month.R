rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )

# Subset salinity and temperature data
saldat <- epcwq3[ which( epcwq3$param=="Sal_top" & year(epcwq3$date) >= 2000 ), ]
tempdat <- epcwq3[ which( epcwq3$param=="Temp_top" & year(epcwq3$date) >= 2000 ), ]
saldat$month <- floor_date( saldat$date, 'month' )
tempdat$month <- floor_date( tempdat$date, 'month' )

# Assemble OTB segment data
saldat.otb <- saldat |> group_by(month) |> dplyr::summarise( sal = mean(value) ) |> as.data.frame()
tempdat.otb <- tempdat |> group_by(month) |> dplyr::summarise( temp = mean(value) ) |> as.data.frame()


# Initialize salinity plot
png( "../figs/temp_sal_otb_month.png", width = 8, height = 8, res = 500, units = 'in' )
par( mfrow = c(2,1), mar = c(3,3.5,1,1) )
plot( sal ~ month, data = saldat.otb,
      ylim = c(0,40), xlab = "", ylab = "", las = 1, 
      col = rgb(0,0,0,0) )
mtext( "Salinity", side = 3, adj = 0, line = 0, cex = 1.2 )
mtext( "ppt", side = 2, line = 2.5, cex = 1 )
abline( h = seq( min(axTicks(2)), max(axTicks(2)),
                 length.out = (length(axTicks(2))-1)*2+1 ),
        col = rgb(0,0,0,0.1) )
abline( v = seq.Date( floor_date(min(saldat.otb$month),'year'),
                      ceiling_date(max(saldat.otb$month),'year'), 'year' ),
                      col = rgb(0,0,0,0.1) )
# Loop through sub-segments to add curves to salinity plot
subsegs <- c("NW","NE","CW","CE","SW","SE")
colors <- hcl.colors( length(subsegs), "mako", alpha = 0.7 )
salcor <- c()
for( i in 1:length(subsegs) ){
  this <- saldat[ which(saldat$subseg==subsegs[i]), ]
  this <- this |> group_by(month) |> dplyr::summarise( sal = mean(value) ) |> as.data.frame()
  lines( sal ~ month, data = this, col = colors[i], lwd = 2 )
  salcordat <- inner_join( this, saldat.otb, by = 'month' )
  salcor[i] <- cor( salcordat$sal.x, salcordat$sal.y )
}  # // end salinity sub-segment loop
legend( 'bottomleft', bty= 'n', ncol = 2,
        legend = paste(subsegs,"sub-segment"),
        lwd = 3, col = colors, seg.len = 1 )



# Initialize temperature plot
plot( temp ~ month, data = tempdat.otb,
      ylim = c(0,40), xlab = "", ylab = "", las = 1, 
      col = rgb(0,0,0,0) )
mtext( "Temperature", side = 3, adj = 0, line = 0, cex = 1.2 )
mtext( "C", side = 2, line = 2.5, cex = 1 )
abline( h = seq( min(axTicks(2)), max(axTicks(2)),
                 length.out = (length(axTicks(2))-1)*2+1 ),
        col = rgb(0,0,0,0.1) )
abline( v = seq.Date( floor_date(min(tempdat.otb$month),'year'),
                      ceiling_date(max(tempdat.otb$month),'year'), 'year' ),
        col = rgb(0,0,0,0.1) )
# Loop through sub-segments to add curves to salinity plot
colors <- hcl.colors( length(subsegs), "redor", alpha = 0.7 )
tempcor <- c()
for( i in 1:length(subsegs) ){
  this <- tempdat[ which(tempdat$subseg==subsegs[i]), ]
  this <- this |> group_by(month) |> dplyr::summarise( temp = mean(value) ) |> as.data.frame()
  lines( temp ~ month, data = this, col = colors[i], lwd = 2 )
  tempcordat <- inner_join( this, tempdat.otb, by = 'month' )
  tempcor[i] <- cor( tempcordat$temp.x, tempcordat$temp.y )
}  # // end salinity sub-segment loop
legend( 'bottomleft', bty= 'n', ncol = 2,
        legend = paste(subsegs,"sub-segment"),
        lwd = 3, col = colors, seg.len = 1 )

dev.off()