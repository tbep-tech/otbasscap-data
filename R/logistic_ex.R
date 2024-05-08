library(tidyverse)

# from docs/eval_paradigm.Rmd
load(file = here::here('data/epcwq2.RData'))
load(file = here::here('data/loads.RData'))

# combine loads and epcwqw
tomod <- epcwq2 |> 
  filter(param %in% c('Chla', 'Turbidity')) |> 
  dplyr::summarise(
    value = mean(value, na.rm = T), 
    .by = c(date, param)
  ) |>
  pivot_wider(names_from = param, values_from = value) |>
  inner_join(loads, by = 'date') |> 
  mutate(chlamet = ifelse(Chla > 9.3, 0, 1)) |> 
  drop_na()

# create glm of prop of chlorophyll exceeding 9.3 given loads and turbidity
mod <- glm(chlamet ~ value + Turbidity, data = tomod, family = 'binomial')

# get model prediction grid and predictions
toprd <- expand_grid(
  value = seq(0, max(tomod$value), length.out = 100),
  Turbidity = c( quantile(tomod$Turbidity, 0.95), mean(tomod$Turbidity) )
)

prds <- predict(mod, type = 'response', newdata = toprd, se.fit = T)
toplo <- toprd |> 
  mutate(
    prd = prds$fit,
    hival = prds$fit + 1.96 * prds$se.fit,
    loval = prds$fit - 1.96 * prds$se.fit, 
    Turbidity = factor( Turbidity, labels = c('95th %tile','mean'), levels = unique(Turbidity) )
  ) 

# get lines to show connection between hypothetical target and certainty of meeting threshold
trgs <- expand_grid(
  value  = c(50), 
  Turbidity = c( quantile(tomod$Turbidity, 0.95), mean(tomod$Turbidity) )
)
lnprds <- predict(mod, type = 'response', newdata = trgs, se.fit = T)
tolns <- trgs |> 
  mutate(
    prd = lnprds$fit,
    hival = lnprds$fit + 1.96 * lnprds$se.fit,
    loval = lnprds$fit - 1.96 * lnprds$se.fit, 
    Turbidity = factor( Turbidity, labels = c('95th %tile','mean'), levels = unique(Turbidity) )
  )

# plot
p <- ggplot(toplo, aes(x = value, y = prd, fill = Turbidity, color = Turbidity, group = Turbidity)) +
  geom_ribbon(aes(ymin = loval, ymax = hival), alpha = 0.2) +
  geom_line() +
  geom_segment(data = tolns, aes(x = value, xend = value, y = 0, yend = hival), linetype = 'dashed', inherit.aes = F) +
  geom_segment(data = tolns, aes(x = 0, xend = value, y = loval, yend = loval), linetype = 'dashed') +  
  geom_segment(data = tolns, aes(x = 0, xend = value, y = hival, yend = hival), linetype = 'dashed') +
  theme_minimal() +
  theme(legend.position = 'top') +
  labs(
    x = 'Monthly TN load', 
    y = 'Probability of attaining 9.3 µg/L',
    caption = 'Hypothetical example of how to use model to set TN load targets'
  ) 

png(here::here('figs/logisticex.png'), width = 6, height = 4, units = 'in', res = 300)
print(p)
dev.off()
