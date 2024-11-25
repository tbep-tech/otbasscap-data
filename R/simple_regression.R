library(dplyr)
library(lubridate)
library(ggplot2)
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
    tn_load_3mo = tn_load + lag(tn_load) + lag(tn_load, 2)
  ) |> 
  filter(date >= as.Date('2000-01-01'))
  

mod <- lm(chla ~ tn_load_3mo, data = chldat)

# TN at 9.3 chla
tn_load_3mo_9.3 <- (9.3 - mod$coefficients[1]) / mod$coefficients[2]

# TN at 8.5 chla
tn_load_3mo_8.5 <- (8.5 - mod$coefficients[1]) / mod$coefficients[2]

p <- ggplot(chldat, aes(x = tn_load_3mo, y = chla)) +
  geom_point(size = 0.5) +
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
    panel.grid.minor = element_blank()
  ) + 
  labs(
    x = 'TN load (3 month rolling sum)',
    y = 'Monthly chlorophyll-a (ug/L)', 
    title = 'OTB data, 2000 to 2021',
    subtitle = paste0('9.3 ug/L chl estimated at ', round(tn_load_3mo_9.3, 1), ' tons TN, ', round(4 * tn_load_3mo_9.3, 1), ' tons TN annual\n8.5 ug/L chl estimated at ', round(tn_load_3mo_8.5, 1), ' tons TN, ', round(4 * tn_load_3mo_8.5, 1), ' tons TN annual'),
    caption = 'Source: EPCHC, Janicki Env.'
  )

png(here::here('figs/simpl_regression.png'), height = 5, width = 5, units = 'in', res = 300)
print(p)
dev.off()
