
rm(list=ls(all=TRUE)) 

library(dplyr)
library(plyr)
library(lubridate)
library(ggplot2)

load('../data-clean/epcwq_clean.RData')
load('../data/loads.RData')

# Specify subsegment
subseg <- unique(epcwq3$subseg)

# Monthly chlorophyll data
chl <- epcwq3[ which( epcwq3$param=="Chla" &
                        epcwq3$subseg %in% subseg ),
               c('date','value') ]
colnames( chl ) <- c('date','chl')
chl <- chl |> group_by(date) |> dplyr::summarise(chl=mean(chl)) |> 
  data.frame()
chl$year <- year( chl$date )

# Annual chlorophyll data 
chl_yr <- chl |> group_by(year) |> dplyr::summarise(chl=mean(chl)) |> 
  as.data.frame()
chl$chl_year_t1 <- mapvalues( chl$year, chl_yr$year, chl_yr$chl )


# Monthly TN loads
loads <- loads[ which( loads$param=="TN load" &
                         year(loads$date) >= 2000 ),
                c('date','value') ]
colnames( loads ) <- c('date','TN_load')
loads$year <- year( loads$date )

# Annual TN loads
loads_yr <- loads |> group_by(year) |> dplyr::summarise(TN_load=sum(TN_load)) |> 
  as.data.frame()
loads$TN_load_year_t1 <- mapvalues( loads$year, loads_yr$year, loads_yr$TN_load )


# Input data
dat <- inner_join( chl[,c('date','chl','chl_year_t1')],
                   loads[,c('date','TN_load','TN_load_year_t1')], 
                   by = 'date' )


# Logistic model
# Specify chlorophyll target
chl_target <- 6.9

# Create binary response variable
dat$y <- ifelse( dat$chl > chl_target, 0, 1 )
# Fit the model
mod <- glm( y ~ TN_load, #+ chl_year_t1 + TN_load_year_t1,
            data = dat, family = 'binomial' )

# Generate predictions
toprd <- expand.grid(
  chl_year_t1 = unique( dat$chl_year_t1 ),
  TN_load = seq( 0, max(dat$TN_load), 0.5 ),
  TN_load_year_t1 = unique( dat$TN_load_year_t1 )
)
prds <- predict( mod, type = 'response',
                 newdata = toprd, se.fit = TRUE )
toplo <- toprd |> 
  mutate(
    prd = prds$fit,
    hival = prds$fit + 1.96 * prds$se.fit,
    loval = prds$fit - 1.96 * prds$se.fit
  )

# Summarise predictions
plot_data <- toplo |> group_by(TN_load) |>
  dplyr::summarise(median=median(prd)) |> as.data.frame()
plot_data$upr <- toplo |> group_by(TN_load) |>
  dplyr::summarise(upr=median(hival)) |> as.data.frame() |> select(upr) 
plot_data$lwr <- toplo |> group_by(TN_load) |>
  dplyr::summarise(lwr=median(loval)) |> as.data.frame() |> select(lwr) 

# Plot predictions
png( "../figs/logistic_model_otb_69.png", width = 7, height = 7, units = 'in', res = 600 )
par(mar=c(4.5,4,2,1))
plot( median ~ TN_load, data = plot_data, type = 'l',
      col = rgb(0.3,0.7,1,1), lwd = 2,
      ylim = c(0,1), xlim = c(0,150), las = 1,
      main = paste0("Probability of attaining ",chl_target," ug/L target (monthly)"),
      ylab = "Probability", yaxt = 'n',
      xlab = "TN load (tons/monthly)", xaxt = 'n' )
axis( 1, at = seq(0,200,20) )
axis( 2, at = seq(0,1,0.2), labels = paste0( seq(0,1,0.2)*100,"%"), las = 1 )
polygon( x = c( plot_data$TN_load, rev(plot_data$TN_load) ),
         y = c( plot_data$upr$upr, rev(plot_data$lwr$lwr) ),
         col = rgb(0.3,0.7,1,0.2), border = rgb(0,0,0,0) )
abline( v = seq(0,200,10), col = rgb(0,0,0,0.1) )  
abline( h = seq(0,1,0.1), col = rgb(0,0,0,0.1) )  
lines( upr$upr ~ TN_load, data = plot_data, col = rgb(0.3,0.7,1,1) )  
lines( lwr$lwr ~ TN_load, data = plot_data, col = rgb(0.3,0.7,1,1) )  
# abline( v = 0 )
# Plot targets
TN_targets <- c(40.5,73)
seg_colors <- c( rgb(0.1,0.6,1,0.7), rgb(0,0.1,0.1,0.7) )
text_pos <- c(-3,3)
for( i in 1:length(TN_targets) ){
  segments( x0 = TN_targets[i], y0 = -1,
            y1 = plot_data$upr$upr[which(plot_data$TN_load==TN_targets[i])],
            lty = 2, col = seg_colors[i] )
  segments( x0 = -100, x1 = TN_targets[i],
            y0 = plot_data$upr$upr[which(plot_data$TN_load==TN_targets[i])],
            lty = 2, col = seg_colors[i]
  )
  segments( x0 = -100, x1 = TN_targets[i],
            y0 = plot_data$median[which(plot_data$TN_load==TN_targets[i])],
            lty = 2, col = seg_colors[i]
  )
  segments( x0 = -100, x1 = TN_targets[i],
            y0 = plot_data$lwr$lwr[which(plot_data$TN_load==TN_targets[i])],
            lty = 2, col = seg_colors[i]
  )
  text( x = rep(text_pos[i],3),
        y = c( plot_data$median[which(plot_data$TN_load==TN_targets[i])],
               plot_data$upr$upr[which(plot_data$TN_load==TN_targets[i])],
               plot_data$lwr$lwr[which(plot_data$TN_load==TN_targets[i])]
        ),
        labels = round( c( 100*plot_data$median[which(plot_data$TN_load==TN_targets[i])],
                           100*plot_data$upr$upr[which(plot_data$TN_load==TN_targets[i])],
                           100*plot_data$lwr$lwr[which(plot_data$TN_load==TN_targets[i])]
        ),0 ), cex = 0.7, offset = 0.1, pos = 3, col = seg_colors[i]
  )
}
legend( 'topright', bty = 'n', lty = 2, col = seg_colors,
        legend = paste0(TN_targets," tons/month")
)
dev.off()
