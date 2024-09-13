rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )
load( "../data/Pyro.Rdata")

# Subset routine Pyro samples
pyro <- pyro[ which(pyro$routine==TRUE), ]

png( "../figs/cumu_pyro-chl_demo.png", width = 10, height = 5, units = 'in', res = 500 )
par(mfrow=c(1,2), mar=c(4,3,3,2))

# Assemble data for analysis
# chlorophyll data (EPC)
epcwq3.sub <- epcwq3[ which( epcwq3$param=="Chla" &
                               year(epcwq3$date) >= 2012 ), ]
epcwq3.sub$month <- floor_date( epcwq3.sub$date, unit = 'month' ) 
chldat <- epcwq3.sub |> group_by(month) |> dplyr::summarise( chl = mean(value) ) |> as.data.frame()
# pyro data (FWC)
pyro.sub <- pyro[ which( pyro$yr>=2012 ), ]
pyro.sub$month <- floor_date( pyro.sub$date, unit = 'month' ) 
pyro.sub <- pyro.sub[ which(complete.cases(pyro.sub)), ]  # pyro==NA means zero cells/L
pyro.sub$logval <- log10( pyro.sub$pyro )
pyrodat <- pyro.sub |> group_by(month) |> dplyr::summarise( pyro = max(logval) ) |> as.data.frame()
# join pyro and chl data by month
pcdat <- inner_join( pyrodat, chldat, by = 'month' )

# Plot histogram of max pyro assoc with mean chl-a at or below 8.5 ug/L
hist( pcdat$pyro[ which(pcdat$chl <= 8.5) ], breaks = 20, freq = TRUE,
      border = rgb(1,1,1,1), col = rgb(0,0,0,0.4),
      main = "", xlab = "", ylab = "",
      las = 1, xaxt = 'n', xlim = c(1,7)
      )
# mtext( expression(atop("(a) Monthly maximum "*italic(P.~bahamense)*" abundance",
#                        "associated with monthly mean chlorophyll-a \u22648.5 ug/L")),
#        side = 3, cex = 1.1 )
mtext( expression("(a) Monthly maximum "*italic(P.~bahamense)*" distribution"),
       side = 3, line = 1.5, adj = 0 )
mtext( expression("      associated with monthly mean chlorophyll-a \u22648.5 ug/L"),
       side = 3, line = 0.5, adj = 0 )
mtext( expression(italic(P.~bahamense)*" (cells/L)"),
       side = 1, line = 2.5 )
mtext( "Number of samples", side = 2, line = 2 )
axis( 1, at = 1:7, labels = c( expression(10^1), expression(10^2),
                               expression(10^3), expression(10^4),
                               expression(10^5), expression(10^6),
                               expression(10^7) ) )
abline( v = pcdat$pyro[ which(pcdat$chl <= 8.5) ] |> median(), 
        col = rgb(1,0.1,0.1,0.6), lwd = 3 )
abline( v = pcdat$pyro[ which(pcdat$chl <= 8.5) ] |> quantile(0.25),
        lty = 3, col = rgb(1,0.1,0.1,0.6), lwd = 3 )
abline( v = pcdat$pyro[ which(pcdat$chl <= 8.5) ] |> quantile(0.75),
        lty = 3, col = rgb(1,0.1,0.1,0.6), lwd = 3 )

# Define function to summarize pyro distribution
distn <- function( x, max_chl ){
  # subset pyro data associated with chl values at or below max_chl
  this <- x$pyro[ which( x$chl <= max_chl ) ]
  # assemble output statistics
  out <- data.frame( median = median(this),
                     min = min(this),
                     lwr_iqr = quantile(this, 0.25),
                     upr_iqr = quantile(this, 0.75),
                     max = max(this)
  )
  return( out )
}  # // end distn()


# Calculate pyro distn statistics
# initiate dataframe
pyro_chl <- data.frame( chl = seq(4,25,0.1),
                        median = NA,
                        min = NA,
                        lwr_iqr = NA,
                        upr_iqr = NA,
                        max = NA
)
# populate rows
for( i in 1:nrow(pyro_chl) ){
  pyro_chl[i,2:6] <- distn( pcdat, max_chl = pyro_chl[i,1] )
}  # // end i loop

# Plot pyro distn as a function of chl_max
plot( chl ~ median, data = pyro_chl, las = 1,
      type = 'l', lwd = 2, col = rgb(0,0.2,0.6,0.7),
      main = "", xlab = "", ylab = "",
      xlim = c(1,7), xaxt = 'n'
      )
mtext( expression("(b) Conditional distribution plot"),
       side = 3, adj = 0, line = 1.5 )
mtext( expression(italic(P.~bahamense)*" (cells/L)"),
       side = 1, line = 2.5 )
mtext( "Chlorophyll a (ug/L)", side = 2, line = 2 )
axis( 1, at = 1:7, las = 1,
      labels = c( expression(10^1), expression(10^2),
                  expression(10^3), expression(10^4),
                  expression(10^5), expression(10^6),
                  expression(10^7) ) )
abline( v = log10( c(1e0,2e0,3e0,4e0,5e0,6e0,7e0,8e0,9e0,
                     1e1,2e1,3e1,4e1,5e1,6e1,7e1,8e1,9e1,
                     1e2,2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,
                     1e3,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,
                     1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,
                     1e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,
                     1e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,
                     1e7,2e7,3e7,4e7,5e7,6e7,7e7,8e7,9e7) ), col = rgb(0,0,0,0.1) )
abline( h = min(pyro_chl$chl):max(pyro_chl$chl), col = rgb(0,0,0,0.1) )
polygon( x = c( pyro_chl$min, rev(pyro_chl$max) ),
         y = c( pyro_chl$chl, rev(pyro_chl$chl) ),
         col = rgb(0,0.2,0.6,0.1), border = rgb(0,0,0,0) )
polygon( x = c( pyro_chl$lwr_iqr, rev(pyro_chl$upr_iqr) ),
         y = c( pyro_chl$chl, rev(pyro_chl$chl) ),
         col = rgb(0,0.2,0.6,0.2), border = rgb(0,0,0,0) )
segments( x0 = 0, x1 = pyro_chl$upr_iqr[which(pyro_chl$chl==8.5)],
          y0 = 8.5,
          lty = 3, col = rgb(1,0.1,0.1,0.6), lwd = 3 )
text( x = 2, y = 9, col = rgb(1,0.1,0.1,0.8), labels = "8.5 ug/L", pos = 3 )
segments( x0 = pyro_chl$median[which(pyro_chl$chl==8.5)],
          y0 = 0, y1 = 8.5,
          lty = 1, col = rgb(1,0.1,0.1,0.6), lwd = 3 )
segments( x0 = pyro_chl$lwr_iqr[which(pyro_chl$chl==8.5)],
          y0 = 0, y1 = 8.5,
          lty = 3, col = rgb(1,0.1,0.1,0.6), lwd = 3 )
segments( x0 = pyro_chl$upr_iqr[which(pyro_chl$chl==8.5)],
          y0 = 0, y1 = 8.5,
          lty = 3, col = rgb(1,0.1,0.1,0.6), lwd = 3 )

segments( x0 = 2.3, x1 = 2.3, y0 = 5.8, y1 = 6.6, lwd = 2 )
arrows( x0 = 2.3, x1 = 2.7, y0 = 5.8, y1 = 5.8, lwd = 2, length = 0.05 )
arrows( x0 = 2.3, x1 = 3.2, y0 = 6.6, y1 = 6.6, lwd = 2, length = 0.05 )
text( x = 1.8, y = 6.2, labels = "target\nrange" )

legend( 'topleft', bty = 'n',
        legend = c("Median",
                   "IQR",
                   "Range"), text.font = 2,
        text.col = c( rgb(0,0.2,0.6,0.9), rgb(0,0.2,0.6,0.5), rgb(0,0.2,0.6,0.3) ) )


dev.off()