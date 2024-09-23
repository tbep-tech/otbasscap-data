rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )

# Specify subsegment and target or threshold value
subsegment <- c("NW","CW")
summer_months <- 6:10  # numeric (e.g. June is 6)
critval <- c(8.5,9.3)  # mgmt target = 8.5 ug/l; NNC threshold = 9.3 ug/l

png( "../figs/chl_attain_summer.png", width = 10, height = 8, units = 'in', res = 600 )
par(mfrow=c(2,2),mar=c(4,5,3,1))

for( i in 1:length(critval) ){
  for( j in 1:length(subsegment) ){
    # Compute annual mean chl conc and attainment
    chl <- epcwq3[ which( epcwq3$param=="Chla" ),
                   c('date','value','subseg') ]
    chl$year <- year( chl$date )
    dat <- chl |> group_by(year) |> dplyr::summarise(value=mean(value)) |>
      as.data.frame()
    colnames( dat ) <- c("year","annual_mean")
    dat$attain <- ifelse( dat$annual_mean < critval[i], 1, 0 )
    
    # Compute summer mean chl conc
    dat <- chl[ which( (month(chl$date) %in% summer_months) & 
                       (chl$subseg %in% subsegment[j]) ), ] |>
      group_by(year) |> dplyr::summarise(summer_mean=mean(value)) |> 
      as.data.frame() |> right_join(dat,'year')
    
    # Fit logistic model to predict probability of attainment from summer means
    mod2 <- glm( attain ~ summer_mean, data = dat, family = 'binomial' )
    prd.x <- seq(0,50,0.1)
    prd.y <- predict( mod2, newdata = data.frame(summer_mean=prd.x),
                      type = 'response', se.fit = TRUE)
    
    # Plot probability of chl attainment
    plot( attain ~ summer_mean, data = dat, yaxt = 'n',
          main = paste0("\nProbability of annual attainment (<",critval[i]," ug/L)"),
          ylab = 'Probability', xlim = c(1,40), ylim = c(0,1),
          xlab = paste0("Summer mean chlorophyll-a (ug/L), ",subsegment[j]," sub-segment"),
          cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3,
          col = rgb(1,0.3,0.2,0.5), pch = 16 )
    if(i==1){
      mtext( paste0("(",letters[j],") ",subsegment[j]," sub-segment"),
             line = 2, side = 3, adj = 0, font = 2 )
    }
    abline( v = seq(0,50,1), col = rgb(0,0,0,0.1) )
    abline( v = seq(0,50,5), col = rgb(0,0,0,0.1) )
    abline( h = seq(0,1,0.1), col = rgb(0,0,0,0.1) )
    axis( 2, at = seq(0,1,0.1), las = 1,
          label = paste0(100*seq(0,1,0.1),"%") )
    lines( prd.y$fit ~ prd.x, col = rgb(0.3,0.4,0.7,1), lwd = 3 )
    polygon( x = c( prd.x, rev(prd.x)), 
             y = c( prd.y$fit+1.96*prd.y$se.fit, rev(prd.y$fit-1.96*prd.y$se.fit) ),
             col = rgb(0.2,0.5,0.9,0.2), border = rgb(0,0,0,0) )
    
  }  # // end subsegment loop j
}  # end critval loop i

dev.off()