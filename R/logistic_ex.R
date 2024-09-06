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
  inner_join(loads[which(loads$param=="TN load"),], by = 'date') |> 
  mutate(chlamet = ifelse(Chla > 9.3, 0, 1)) |> 
  drop_na()

# create glm of prop of chlorophyll exceeding 9.3 given loads and turbidity
mod <- glm(chlamet ~ value + Turbidity, data = tomod, family = 'binomial')

# get model prediction grid and predictions
toprd <- expand_grid(
  value = seq(0, max(tomod$value), length.out = 100),
  Turbidity = c( quantile(tomod$Turbidity, 0.99), mean(tomod$Turbidity) )
)

prds <- predict(mod, type = 'response', newdata = toprd, se.fit = T)
toplo <- toprd |> 
  mutate(
    prd = prds$fit,
    hival = prds$fit + 1.96 * prds$se.fit,
    loval = prds$fit - 1.96 * prds$se.fit, 
    Turbidity = factor( Turbidity, labels = c('Elevated','Normal'), levels = unique(Turbidity) )
  ) 

# get lines to show connection between hypothetical target and certainty of meeting threshold
trgs <- expand_grid(
  value  = c(50), 
  Turbidity = c( quantile(tomod$Turbidity, 0.99), mean(tomod$Turbidity) )
)
lnprds <- predict(mod, type = 'response', newdata = trgs, se.fit = T)
tolns <- trgs |> 
  mutate(
    prd = lnprds$fit,
    hival = lnprds$fit + 1.96 * lnprds$se.fit,
    loval = lnprds$fit - 1.96 * lnprds$se.fit, 
    Turbidity = factor( Turbidity, labels = c('Elevated','Normal'), levels = unique(Turbidity) )
  )

# plot
p <- ggplot(toplo, aes(x = value, y = prd, group = Turbidity)) +
  geom_ribbon(aes(ymin = loval, ymax = hival, fill = Turbidity), alpha = 0.2) +
  geom_line(aes(color = Turbidity)) +
  geom_segment(data = tolns, aes(x = value, xend = value, y = 0, yend = prd), linetype = 'dashed', inherit.aes = F) +
  geom_segment(data = tolns, aes(color = Turbidity, x = 0, xend = value, y = prd, yend = prd), linetype = 'dashed') +  
  # geom_segment(data = tolns, aes(x = 0, xend = value, y = hival, yend = hival), linetype = 'dashed') +
  scale_color_manual(values = c('tomato1', 'dodgerblue')) +
  scale_fill_manual(values = c('tomato1', 'dodgerblue')) +
  coord_cartesian(xlim = c(0, 200), ylim = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    legend.position = 'top',
    panel.background = element_rect(fill = 'transparent', color = NA),
    plot.background = element_rect(fill = 'transparent', color = NA),
    legend.background = element_rect(fill = 'transparent', color = NA),
    legend.box.background = element_rect(fill = 'transparent', color = NA)
  ) +
  labs(
    x = 'Monthly nitrogen load', 
    fill = 'Climate condition', 
    color = 'Climate condition',
    y = 'Probability of meeting water quality goal'
  )

png(here::here('figs/logisticex.png'), width = 4.5, height = 4, units = 'in', res = 300, bg = 'transparent')
print(p)
dev.off()
