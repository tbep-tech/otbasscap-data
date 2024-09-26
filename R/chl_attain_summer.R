rm(list=ls(all=TRUE)) 

library(dplyr)
library(lubridate)

# Load data from file
load( "../data-clean/epcwq_clean.RData" )

# Specify subsegment and target or threshold value
subsegment <- c("NW","CW")
summer_months <- 6:10  # numeric (e.g. June is 6)
critval <- c(8.5,9.3)  # mgmt target = 8.5 ug/l; NNC threshold = 9.3 ug/l

# Initialize tables to hold summertime concs assoc with 50% and 75% prob
prob_50 <- matrix(NA,nrow=2,ncol=2)
  row.names(prob_50) <- rep(NA,2)
  colnames(prob_50) <- rep(NA,2)
prob_75 <- matrix(NA,nrow=2,ncol=2)
  row.names(prob_75) <- rep(NA,2)
  colnames(prob_75) <- rep(NA,2)

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
    
    # Print to predicted values to console
    # prob <- 0.50
    # df <- data.frame( x = prd.x,
    #                   y_low = prd.y$fit-1.96*prd.y$se.fit,
    #                   y_med = prd.y$fit,
    #                   y_hi  = prd.y$fit+1.96*prd.y$se.fit )
    # out <- df[ which( (df$y_low>(prob*0.95) & df$y_low<(prob*1.05)) |
    #                   (df$y_med>(prob*0.95) & df$y_med<(prob*1.05)) |
    #                   (df$y_hi>(prob*0.95) & df$y_hi<(prob*1.05)) ), ]
    # cat( "\n",subsegment[j],"subseg    critval:",critval[i],"   prob:",prob,"\n" ); print( out )
    
    # Populate 50% and 75% probability tables
    df_50 <- data.frame( x = prd.x,
                         y_low = prd.y$fit-1.96*prd.y$se.fit,
                         y_med = prd.y$fit,
                         y_hi  = prd.y$fit+1.96*prd.y$se.fit )
    prob_50[i,j] <- paste( format(df_50$x[ which( df_50$y_low<0.50 )[1] ],nsmall=1),
                           format(df_50$x[ which( df_50$y_med<0.50 )[1] ],nsmall=1),
                           format(df_50$x[ which( df_50$y_hi<0.50 )[1] ],nsmall=1),
                           sep = " / " )
    row.names(prob_50)[i] <- paste0(critval[i]," ug/L")
    colnames(prob_50)[j] <- paste0(subsegment[j]," sub-segment")
    
    df_75 <- data.frame( x = prd.x,
                         y_low = prd.y$fit-1.96*prd.y$se.fit,
                         y_med = prd.y$fit,
                         y_hi  = prd.y$fit+1.96*prd.y$se.fit )
    prob_75[i,j] <- paste( format(df_75$x[ which( df_50$y_low<0.75 )[1] ],nsmall=1),
                           format(df_75$x[ which( df_50$y_med<0.75 )[1] ],nsmall=1),
                           format(df_75$x[ which( df_50$y_hi<0.75 )[1] ],nsmall=1),
                           sep = " / " )
    row.names(prob_75)[i] <- paste0(critval[i]," ug/L")
    colnames(prob_75)[j] <- paste0(subsegment[j]," sub-segment")
    
    # Add lines to plots to indicate 75% probability and associated concentration
    if( i == 1 ){
    segments( x0 = -20, x1 = df_75$x[ which( df_50$y_low<0.75 )[1] ], y0 = 0.75, lty = 2 )
    segments( x0 = df_75$x[ which( df_50$y_low<0.75 )[1] ], y0 = -1, y1 = 0.75, lty = 2 )
    }
    
  }  # // end subsegment loop j
}  # end critval loop i

dev.off()

# Export probability tables
write.csv( prob_50, "../data/chl_attain_summer_50.csv", row.names = TRUE )
write.csv( prob_75, "../data/chl_attain_summer_75.csv", row.names = TRUE )