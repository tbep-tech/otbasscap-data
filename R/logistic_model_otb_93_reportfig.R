
png( "../figs/logistic_model_otb_93_reportfig.png", width = 5, height = 10, units = 'in', res = 600 )
par(mar=c(4.5,4,2,1),mfrow=c(2,1))

# logistic model for absolute TMDL ----------------------------------------


rm(list=ls(all=TRUE)) 

library(dplyr)
library(plyr)
library(lubridate)
library(ggplot2)

load(('../data-clean/epcwq_clean.RData'))
load(('../data/loads.RData'))

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
chl_target <- 9.3

# Create binary response variable
dat$y <- ifelse( dat$chl > chl_target, 0, 1 )
# Fit the model
mod <- glm( y ~ TN_load, #+ chl_year_t1 + TN_load_year_t1,
            data = dat, family = 'binomial' )

# Generate predictions
toprd <- data.frame( TN_load = seq( 0, max(dat$TN_load), 0.5 ) )
# toprd <- expand.grid(
#   chl_year_t1 = unique( dat$chl_year_t1 ),
#   TN_load = seq( 0, max(dat$TN_load), 0.5 ),
#   TN_load_year_t1 = unique( dat$TN_load_year_t1 )
#   )
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
plot( median ~ TN_load, data = plot_data, type = 'l',
      col = rgb(0.3,0.7,1,1), lwd = 2,
      ylim = c(0,1), xlim = c(0,150), las = 1,
      main = "",
      ylab = "Probability", yaxt = 'n',
      xlab = "", xaxt = 'n' )
mtext( "(a) Chl-a threshold attainment and TN load",
       side = 3, line = 0, cex = 1.2, adj = 0 )
mtext( "TN load (tons/month)", side = 1, line = 2, cex = 1 )
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
            lty = 2, col = seg_colors[i], lwd = 1.5
            )
  segments( x0 = -100, x1 = TN_targets[i],
            y0 = plot_data$upr$upr[which(plot_data$TN_load==TN_targets[i])],
            lty = 2, col = seg_colors[i], lwd = 1.5
  )
  segments( x0 = -100, x1 = TN_targets[i],
            y0 = plot_data$median[which(plot_data$TN_load==TN_targets[i])],
            lty = 2, col = seg_colors[i], lwd = 1.5
  )
  segments( x0 = -100, x1 = TN_targets[i],
            y0 = plot_data$lwr$lwr[which(plot_data$TN_load==TN_targets[i])],
            lty = 2, col = seg_colors[i], lwd = 1.5
  )
  # text( x = rep(text_pos[i],3),
  #       y = c( plot_data$median[which(plot_data$TN_load==TN_targets[i])],
  #              plot_data$upr$upr[which(plot_data$TN_load==TN_targets[i])],
  #              plot_data$lwr$lwr[which(plot_data$TN_load==TN_targets[i])]
  #              ),
  #       labels = round( c( 100*plot_data$median[which(plot_data$TN_load==TN_targets[i])],
  #                          100*plot_data$upr$upr[which(plot_data$TN_load==TN_targets[i])],
  #                          100*plot_data$lwr$lwr[which(plot_data$TN_load==TN_targets[i])]
  #                          ),0 ), cex = 0.7, offset = 0.1, pos = 3, col = seg_colors[i]
  #      )
}
legend( 'topright', bty = 'n', lty = 2, lwd = 1.5, col = seg_colors,
        legend = paste0(TN_targets," tons/month")
)



# logistic model for delivery ratio ---------------------------------------


rm(list=ls(all=TRUE))

load(('../data-clean/epcwq_clean.RData'))  # WQ data
load(('../data/loads.RData'))  # monthly TN loads
hydro <- read.csv(("../data-raw/TB_hydro_monthly.csv"))  # monthly hydro loads
load(("../data-raw/totanndat.RData"))  # annual TN and hydro loads

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

# Annual mean chlorophyll
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
loads_yr <- totanndat[ which(totanndat$year>=2000 &
                               totanndat$bay_segment=="Old Tampa Bay"),
                       c('year','tn_load')] |> as.data.frame()
loads$TN_load_year_t1 <- mapvalues( loads$year, loads_yr$year, loads_yr$tn_load )


# Monthly hydro loads
hydro <- hydro[ which( hydro$bay_segment=="Old Tampa Bay" &
                         hydro$year >= 2000 ),
                c('year','month','hy_load_106_m3_mo') ]
hydro$month <- paste0( hydro$year,'-',hydro$month,'-01' ) |> as.Date()
colnames( hydro ) <- c('year','date','hydro_load')

# Annual hydro loads
hydro_yr <- totanndat[ which(totanndat$year>=2000 &
                               totanndat$bay_segment=="Old Tampa Bay"),
                       c('year','hy_load')] |> as.data.frame()
# hydro_yr <- hydro |> group_by(year) |> dplyr::summarise(hydro_load=sum(hydro_load)) |>
#               as.data.frame()
hydro$hydro_load_year_t1 <- mapvalues( hydro$year, hydro_yr$year, hydro_yr$hy_load )


# TN delivery ratio
deliv <- inner_join( loads[,c('date','TN_load','TN_load_year_t1')],
                     hydro[,c('date','hydro_load','hydro_load_year_t1')],
                     by = 'date'
)
deliv$ratio <- deliv$TN_load / deliv$hydro_load
deliv$ratio_year_t1 <- deliv$TN_load_year_t1 / deliv$hydro_load_year_t1
deliv <- deliv[ ,c('date','ratio','ratio_year_t1') ]


# Input data
dat <- inner_join( chl[,c('date','chl','chl_year_t1')],
                   deliv[,c('date','ratio','ratio_year_t1')],
                   by = 'date' )


# Logistic model
# Specify chlorophyll target
chl_target <- 9.3

# Create binary response variable
dat$y <- ifelse( dat$chl > chl_target, 0, 1 )
# Fit the model
mod <- glm( y ~ ratio, #+ chl_year_t1 + ratio_year_t1,
            data = dat, family = 'binomial' )

# Generate predictions
toprd <- data.frame( ratio = seq( 0, max(dat$ratio), 0.01 ) )
# toprd <- expand.grid(
#   chl_year_t1 = unique( dat$chl_year_t1 ),
#   ratio = seq( 0, max(dat$ratio), 0.01 ),
#   ratio_year_t1 = unique( dat$ratio_year_t1 )
#   )
prds <- predict( mod, type = 'response',
                 newdata = toprd, se.fit = TRUE )
toplo <- toprd |>
  mutate(
    prd = prds$fit,
    hival = prds$fit + 1.96 * prds$se.fit,
    loval = prds$fit - 1.96 * prds$se.fit
  )

# Summarise predictions
plot_data <- toplo |> group_by(ratio) |>
  dplyr::summarise(median=median(prd)) |> as.data.frame()
plot_data$upr <- toplo |> group_by(ratio) |>
  dplyr::summarise(upr=median(hival)) |> as.data.frame() |> select(upr)
plot_data$lwr <- toplo |> group_by(ratio) |>
  dplyr::summarise(lwr=median(loval)) |> as.data.frame() |> select(lwr)

# Plot predictions
plot( median ~ ratio, data = plot_data, type = 'l',
      col = rgb(0.3,0.7,1,1), lwd = 2,
      ylim = c(0,1), xlim = c(0,1.7), las = 1,
      main = "",
      ylab = "Probability", yaxt = 'n',
      xlab = "", xaxt = 'n' )
mtext( "(b) Chl-a threshold attainment and\nTN delivery ratio",
       side = 3, line = 0, cex = 1.2, adj = 0 )
mtext( "TN delivery ratio (tons/Mm3 per year)", side = 1, line = 2, cex = 1 )
axis( 1, at = seq(0,2,0.2) )
axis( 2, at = seq(0,1,0.2), labels = paste0( seq(0,1,0.2)*100,"%"), las = 1 )
polygon( x = c( plot_data$ratio, rev(plot_data$ratio) ),
         y = c( plot_data$upr$upr, rev(plot_data$lwr$lwr) ),
         col = rgb(0.3,0.7,1,0.2), border = rgb(0,0,0,0) )
abline( v = seq(0,2,0.1), col = rgb(0,0,0,0.1) )
abline( h = seq(0,1,0.1), col = rgb(0,0,0,0.1) )
lines( upr$upr ~ ratio, data = plot_data, col = rgb(0.3,0.7,1,1) )
lines( lwr$lwr ~ ratio, data = plot_data, col = rgb(0.3,0.7,1,1) )
# abline( v = 0 )
# Plot targets
deliv_targets <- c(1.08,1.71)
seg_colors <- c( rgb(0.1,0.6,1,0.7), rgb(0,0.1,0.1,0.7) )
text_pos <- c(-0.03,0.03)
for( i in 1:length(deliv_targets) ){
  segments( x0 = deliv_targets[i], y0 = -1,
            y1 = plot_data$upr$upr[which(plot_data$ratio==deliv_targets[i])],
            lty = 2, col = seg_colors[i], lwd = 1.5
            )
  segments( x0 = -100, x1 = deliv_targets[i],
            y0 = plot_data$upr$upr[which(plot_data$ratio==deliv_targets[i])],
            lty = 2, col = seg_colors[i], lwd = 1.5
  )
  segments( x0 = -100, x1 = deliv_targets[i],
            y0 = plot_data$median[which(plot_data$ratio==deliv_targets[i])],
            lty = 2, col = seg_colors[i], lwd = 1.5
  )
  segments( x0 = -100, x1 = deliv_targets[i],
            y0 = plot_data$lwr$lwr[which(plot_data$ratio==deliv_targets[i])],
            lty = 2, col = seg_colors[i], lwd = 1.5
  )
  # text( x = rep(text_pos[i],3),
  #       y = c( plot_data$median[which(plot_data$ratio==deliv_targets[i])],
  #              plot_data$upr$upr[which(plot_data$ratio==deliv_targets[i])],
  #              plot_data$lwr$lwr[which(plot_data$ratio==deliv_targets[i])]
  #              ),
  #       labels = round( c( 100*plot_data$median[which(plot_data$ratio==deliv_targets[i])],
  #                          100*plot_data$upr$upr[which(plot_data$ratio==deliv_targets[i])],
  #                          100*plot_data$lwr$lwr[which(plot_data$ratio==deliv_targets[i])]
  #                          ),0 ), cex = 0.7, offset = 0.1, pos = 3, col = seg_colors[i]
  #      )
}
legend( 'topright', bty = 'n', lty = 2, lwd = 1.5, col = seg_colors,
        legend = paste0(deliv_targets," tons/Mm3")
)

dev.off()