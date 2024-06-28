# This script computes correlations between Pyro abundance and chlorophyll-a
# concentrations in the EPC data, for OTB sub-segments.
# 2024 Miles Medina (ECCO Scientific)

rm(list=ls(all=TRUE)) 

library(plyr)
library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )
load( "../data-clean/epcphyto.Rdata")

# Specify subsegment
subsegs <- c("NW","NE","CW","CE","SW","SE")

# Assemble data for plots
# Initalize lists
pcdat <- list()  # pyro and chl values
lim <- list( x = c(NA,NA), y = c(NA,NA) )  # axis limits for plotting
# Populate objects
for( i in 1:length(subsegs) ){
  
  # Assemble data for analysis
  # chlorophyll data (EPC)
  epcwq3.sub <- epcwq3[ which( epcwq3$param=="Chla" &
                                 year(epcwq3$date) >= 2012 &
                                 epcwq3$subseg==subsegs[i] ), ]
  chldat <- epcwq3.sub |> group_by(date) |> summarise( chl = mean(value) ) |> as.data.frame()
  # pyro data (EPC)
  pyro.sub <- phyto[ which( phyto$name == "Pyrodinium bahamense" &
                            phyto$subsegment==subsegs[i] &
                            year(phyto$date)>=2012 ), ]

  pyro.sub <- pyro.sub[ ,c('date','cellcount') ]
  pyro.sub <- pyro.sub[ which(complete.cases(pyro.sub)), ]  # pyro==NA means zero cells/L
  pyro.sub$logval <- log10( pyro.sub$cellcount )
  pyrodat <- pyro.sub |> group_by(date) |> summarise( pyro = max(logval) ) |> as.data.frame()
  # join pyro and chl data by month
  pcdat[[i]] <- inner_join( pyrodat, chldat, by = 'date' )
  # assign name to list item
  names(pcdat)[i] <- subsegs[i]
  # Compute xlim and ylim for plotting
  lim$x <- range( c( lim$x, pcdat[[i]]$chl ), na.rm = TRUE )
  lim$y <- range( c( lim$y, min(pcdat[[i]]$pyro)-1,
                     max(pcdat[[i]]$pyro)+1 ), na.rm = TRUE )
  
}  # // end i loop


# Plotting parameters
png( "../figs/pyro_chl_corr_epc.png", height = 9, width = 9, units = 'in', res = 500 )
par( mfrow=c(3,2), mar=c(5,4,2,1) )

for( i in 1:length(pcdat) ){
  
  # Scatterplot
  plot( pcdat[[i]]$pyro ~ pcdat[[i]]$chl,
        xlim = lim$x, ylim = lim$y,
        main = paste0( names(pcdat)[i]," OTB: Pyrodinium vs. chlorophyll"),
        xlab = "Chlorophyll-a (ug/L)",
        ylab = "Pyro abundance (cells/L)", yaxt = 'n',
        col = rgb(1,0.1,0.2,0.6), lwd = 2,
        cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3 )
  axis( 2, at = 1:7, las = 1, cex.axis = 1.3,
        labels = c( expression(10^1), expression(10^2),
                    expression(10^3), expression(10^4),
                    expression(10^5), expression(10^6),
                    expression(10^7)) )
  abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  
  # Linear model
  mod1 <- lm( pcdat[[i]]$pyro ~ pcdat[[i]]$chl )
  intc <- summary(mod1)$coefficients[1]
  slop <- summary(mod1)$coefficients[2]
  R2 <- summary(mod1)$r.squared |> round(3)
  pval <- summary(mod1)$coefficients[8]
  if( pval > 0.05 ){
    pval.string <- "P > 0.05"
  } else if( pval < 0.001 ){
    pval.string <- "P < 0.001"
  } else if( pval < 0.01 ){
    pval.string <- "P < 0.01"
  } else {
    pval.string <- "P < 0.05"
  }
  
  # Add lm results to plot
  segments( x0 = min(pcdat[[i]]$chl), x1 = max(pcdat[[i]]$chl),
            y0 = intc + slop * min(pcdat[[i]]$chl),
            y1 = intc + slop * max(pcdat[[i]]$chl),
            lty = 2
  )
  legend( "bottom", bty = 'n',
          legend = paste0( "N = ", nrow(pcdat[[i]]), "     ",
                           "R2 = ", R2, "     ", pval.string ),
          cex = 1.3 )
  
}  # // end subseg loop (plotting)

dev.off()