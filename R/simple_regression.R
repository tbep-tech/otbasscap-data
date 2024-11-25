library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(patchwork)
library(tbeptools)

load(here::here('data/loads.RData'))

lddat <- loads |> 
  filter(param == 'TN load') |> 
  select(date, tn_load = value)

chldat <- anlz_avedat(epcdata)$mos |> 
  filter(bay_segment == 'OTB' & var == 'mean_chla') |> 
  mutate(
    date = make_date(yr, mo, 1)
  ) |> 
  select(date, chla = val) |> 
  inner_join(lddat, by = 'date') |> 
  mutate(
    tn_load_3mo = tn_load + lag(tn_load) + lag(tn_load, 2),
    qrt = factor(lubridate::quarter(date, with_year = F), levels = 1:4, labels = c('JFM', 'AMJ', 'JAS', 'OND'))
  ) |> 
  filter(date >= as.Date('2000-01-01'))
  

mod <- lm(chla ~ tn_load_3mo, data = chldat)

# TN at 9.3 chla
tn_load_3mo_9.3 <- (9.3 - mod$coefficients[1]) / mod$coefficients[2]

# TN at 8.5 chla
tn_load_3mo_8.5 <- (8.5 - mod$coefficients[1]) / mod$coefficients[2]

p1 <- ggplot(chldat, aes(x = tn_load_3mo, y = chla)) +
  geom_point(size = 2, aes(color = qrt), alpha = 0.8) +
  geom_smooth(method = 'lm', se = T, formula = y ~ x, color = 'black') +
  geom_segment(
    x = tn_load_3mo_9.3, xend = tn_load_3mo_9.3, y = 0, yend = 9.3,
    linetype = 'dashed'
  ) +
  geom_segment(
    x = tn_load_3mo_8.5, xend = tn_load_3mo_8.5, y = 0, yend = 8.5,
    linetype = 'dashed'
  ) +
  geom_segment(
    x = 0, xend = tn_load_3mo_9.3, y = 9.3, yend = 9.3,
    linetype = 'dashed'
  ) +
  geom_segment(
    x = 0, xend = tn_load_3mo_8.5, y = 8.5, yend = 8.5,
    linetype = 'dashed'
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(), 
    legend.position = 'bottom'
  ) + 
  labs(
    x = 'TN load (3 month rolling sum)',
    y = 'Monthly chlorophyll-a (ug/L)', 
    title = 'OTB data, 2000 to 2021',
    color = 'Season',
    subtitle = paste0('9.3 ug/L chl estimated at ', round(tn_load_3mo_9.3, 1), ' tons TN, ', round(4 * tn_load_3mo_9.3, 1), ' tons TN annual\n8.5 ug/L chl estimated at ', round(tn_load_3mo_8.5, 1), ' tons TN, ', round(4 * tn_load_3mo_8.5, 1), ' tons TN annual'),
    caption = 'Source: EPCHC, Janicki Env.'
  )

stochastic_backcalc <- function(model, y_value, n_draws = 1000) {
  
  # Extract coefficient estimates and standard errors
  coef_estimates <- summary(model)$coefficients
  intercept <- coef_estimates["(Intercept)", "Estimate"]
  intercept_se <- coef_estimates["(Intercept)", "Std. Error"]
  slope <- coef_estimates["tn_load_3mo", "Estimate"]
  slope_se <- coef_estimates["tn_load_3mo", "Std. Error"]
  
  # Generate random draws for intercept and slope
  intercept_draws <- rnorm(n_draws, intercept, intercept_se)
  slope_draws <- rnorm(n_draws, slope, slope_se)
  
  # Calculate x-values for each draw
  x_draws <- (y_value - intercept_draws) / slope_draws
  
  #
  return(x_draws)
  # return(list(
  #   mean = mean(x_draws),
  #   median = median(x_draws),
  #   lower_ci = quantile(x_draws, 0.025),
  #   upper_ci = quantile(x_draws, 0.975),
  #   full_draws = x_draws
  # ))

}

stoc93 <- 
stoc85 <- 

toplo <- tibble(
  stoc9.3 = stochastic_backcalc(mod, 9.3),
  stoc8.5 = stochastic_backcalc(mod, 8.5)
  ) |> 
  pivot_longer(everything(), names_to = 'chla', values_to = 'tn_load_3mo') |> 
  mutate(
    chla = gsub('stoc', '', chla) |> factor(), 
    tn_load = 4 * tn_load_3mo
  )
medv <- toplo |> 
  summarise(
    tn_load = median(tn_load), 
    .by = chla
  )

# create ggplot with points jitter and a horizontal line for the median
p2 <- ggplot(toplo, aes(x = chla, y = tn_load)) + 
  geom_jitter(width = 0.2, height = 0, size = 0.5) +
  geom_boxplot(data = medv, width = 0.4, color = 'red') + 
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  ) + 
  labs(
    x = 'Chlorophyll-a (ug/L)',
    y = 'Annual TN load (tons)',
    title = 'Stochastic back-calculation of TN load\nfor 9.3 and 8.5 ug/L chl',
  ) 

p <- p1 + p2 +  plot_layout(ncol = 2)

png(here::here('figs/simple_regression.png'), height = 5.5, width = 10, units = 'in', res = 300)
print(p)
dev.off()
