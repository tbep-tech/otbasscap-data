rm(list=ls(all=TRUE)) 

library(tidyverse)

# define logistic modeling function for P( y < ytarget ) ~ x1 + x2.
#   Arguments:
#     wqdat, loaddat    dataframes containing water quality and loading data, incl.
#                         columns date, param, value, unit, site
#     y                 response variable (e.g. Chla), as appears in wqdat$param
#     x1                first predictor (e.g. TN load), as appears in loaddat$param
#     x2                second predictor (e.g. temperature), as appears in wqdat$param
#     ytarget           numeric target for response variable y (e.g. 9.3 mg/l Chla)
#     xtarget           numeric target(s) for predictor x1 (e.g. 100 tons TN load)
#     qnt               quantile(s) to plot for predictor x2 (numeric values in [0,1])
#     wqloc             subsegment(s) for subsetting wqdat (default uses all unique wqdat$subsegs)
#     plot              logical, displays plot if TRUE
logim <- function( wqdat, loaddat,
                   y, x1, x2,
                   ytarget, xtarget,
                   qnt, wqloc = unique(wqdat$subseg),
                   plot = TRUE ){
  
  # join wq and load data
  wqsub <- wqdat[ which( wqdat$subseg %in% wqloc ), ]
  wqsub <- wqsub |> 
    filter(param %in% c(y, x2)) |> 
    dplyr::summarise(
      value = mean(value, na.rm = TRUE ), 
      .by = c(date, param) ) |>
    pivot_wider(names_from = param, values_from = value)
  wqsub$date <- wqsub$date |> floor_date('month')
  
  loadsub <- loaddat |>
    filter( param==x1 ) |>
    select( c('date','param','value') ) |> 
    pivot_wider( names_from = param, values_from = value )
  
  tomod <- inner_join( wqsub, loadsub, by = 'date' ) |>
    setNames( c('date','y','x2','x1') ) |>
    mutate( bin = ifelse( y > ytarget, 0, 1 ) ) |> 
    drop_na()
  
  y.unit <- wqdat$unit[ which(wqdat$param==y) ] |> unique()
  
  # logistic model for probability of attaining ytarget
  mod <- glm( bin ~ x1 + x2,
              data = tomod, family = 'binomial' )
  
  # get model prediction grid and predictions
  toprd <- expand_grid(
    x1 = seq( 0, max(tomod$x1), length.out = 100 ),
    x2 = c( quantile( tomod$x2, qnt ), mean(tomod$x2) )
  )
  prds <- predict( mod, type = 'response',
                   newdata = toprd, se.fit = TRUE )
  toplo <- toprd |> 
    mutate(
      prd = prds$fit,
      hival = prds$fit + 1.96 * prds$se.fit,
      loval = prds$fit - 1.96 * prds$se.fit,
      x2 = factor( x2, labels = c( paste0(qnt*100,'th %ile'), 'mean' ),
                   levels = unique(x2) ) 
    )
  
  # get lines to highlight probabilities at specified x1 value
  trgs <- expand_grid(
    x1  = xtarget, 
    x2 = c( quantile(tomod$x2, qnt), mean(tomod$x2) )
  )
  lnprds <- predict(mod, type = 'response',
                    newdata = trgs, se.fit = T)
  tolns <- trgs |> 
    mutate(
      prd = lnprds$fit,
      hival = lnprds$fit + 1.96 * lnprds$se.fit,
      loval = lnprds$fit - 1.96 * lnprds$se.fit, 
      x2 = factor( x2, labels = c( paste0(qnt*100,'th %ile'), 'mean' ),
                   levels = unique(x2) )
    )
  
  # plot
  p <- ggplot( toplo, aes( x = x1,
                           y = prd,
                           fill = x2,
                           color = x2,
                           group = x2)) +
    geom_ribbon( aes( ymin = loval, ymax = hival), alpha = 0.2 ) +
    geom_line() +
    geom_segment( data = tolns, aes(x = x1, xend = x1,
                                    y = 0, yend = hival),
                  linetype = 'dashed', inherit.aes = F ) +
    geom_segment( data = tolns, aes( x = 0, xend = x1,
                                     y = loval, yend = loval),
                  linetype = 'dashed') +  
    geom_segment( data = tolns, aes( x = 0, xend = x1,
                                     y = hival, yend = hival),
                  linetype = 'dashed') +
    theme_minimal() +
    theme( legend.position = 'top' ) +
    labs( x = 'Monthly TN load (tons)', 
          y = paste0('Probability of attaining ',y,' target (',ytarget,' ',y.unit,')') ) +
    coord_cartesian(ylim = c(0,1) ) 
  
  p$labels$fill <- x2
  p$labels$colour <- x2
  
  # display plot
  if( plot ){
    print( p )
  }
  
  # assemble function arguments, data, and predicted values
  inputs <- formals()
  inputs$ytarget <- ytarget
  inputs$xtarget <- xtarget
  inputs$y <- y
  inputs$x1 <- x1
  inputs$x2 <- x2
  inputs$wqloc <- wqloc
  inputs$qnt <- qnt
  
  return( list( plot = p,
                model = mod,
                probs = tolns,
                args = inputs,
                dat = tomod,
                pred = toplo
  )
  )
  
}  # // end logim()


# load data
load('../data-clean/epcwq_clean.RData')
load('../data/loads.RData')

# run models
mod.turb <- logim( wqdat = epcwq3, loaddat = loads,
                   y = "Chla", x1 = "TN load", x2 = "Turbidity",
                   ytarget = 9.3, xtarget = c(40),
                   qnt = c(0.95)
                   # , wqloc = c("NW")
                   )

mod.sal  <- logim( wqdat = epcwq3, loaddat = loads,
                   y = "Chla", x1 = "TN load", x2 = "Sal_top",
                   ytarget = 9.3, xtarget = c(40),
                   qnt = c(0.25,0.75)
                   , wqloc = c("NW")
                   )

mod.temp <- logim( wqdat = epcwq3, loaddat = loads,
                   y = "Chla", x1 = "TN load", x2 = "Temp_top",
                   ytarget = 9.3, xtarget = c(40),
                   qnt = c(0.95)
                   , wqloc = c("SW")
                   )
