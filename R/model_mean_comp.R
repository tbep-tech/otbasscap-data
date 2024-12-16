library(tidyverse)
library(patchwork)

dat <- tibble(
  obser = c(2.32, 2.84, 10.76, 9.18, 5.92, 3.83, 15.03, 8.68, 6.91, 12.50),
  reduc = c(2.32, 2.84, 10.56, 3.69, 7.82, 3.79, 12.57, 6.82, 6.94, 12.39)
)

# summarise means and 95% ci
toplo <- dat |> 
  pivot_longer(cols = c(obser, reduc), names_to = "group") |> 
  mutate(
    group = factor(group, levels = c('obser', 'reduc'), labels = c('Observed tons/yr', '486 tons/yr'))
  ) |> 
  summarise(
    mean = mean(value),
    lower = mean - 1.96 * sd(value) / sqrt(n()),
    upper = mean + 1.96 * sd(value) / sqrt(n()),
    .by = group
  )

p1 <- ggplot(toplo, aes(x = group, y = mean, ymin = lower, ymax = upper)) +
  geom_col() +
  geom_errorbar(width = 0.2) +
  theme_minimal() + 
  theme(
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 12)) +
  labs(
    x = NULL,
    y = "ug/L +/- 95% CI", 
    title = "Mean modeled chl-a 2012-2021",
    subtitle = 'Observed vs TN load reduction'
  )

# run t.test and plot marginal effect
p2 <- t.test(dat$reduc, dat$obser, paired = TRUE) |> 
  broom::tidy() |> 
  ggplot(aes(y = estimate, x = 'Diff')) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() + 
  theme(
    panel.grid.minor = element_blank(), 
    axis.text.x = element_blank(), 
    panel.grid.major.x = element_blank()
  ) +
  labs(
    y = "ug/L +/- 95% CI",
    x = NULL, 
    subtitle = "Difference in means"
  )

p <- p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))

# png('~/Desktop/scenariocomp.png', height = 4, width = 6, res = 300, units = 'in')
# print(p)
# dev.off()
