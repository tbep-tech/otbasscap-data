# setup ---------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(here)

load(here('data-clean/epcwq_clean.RData'))
load(here('data/loads.RData'))
hydro <- read.csv(here("data-raw/TB_hydro_monthly.csv"))

# hydro loads
hyd <- hydro %>%
  filter(bay_segment == "Old Tampa Bay", year >= 2000) %>%
  mutate(date = as.Date(paste0(year, "-", month, "-01"))) %>%
  select(year, date, hydro_load = hy_load_106_m3_mo) %>%
  mutate(
    hydroann_load = sum(hydro_load), 
    .by = year
  )

# prep loading data
lds <- loads %>%
  filter(param == "TN load", year(date) >= 2000) %>%
  select(date, TN_load = value) %>%
  inner_join(hyd, by = 'date') |> 
  mutate(
    ratio = TN_load / hydro_load,
    TN_load_norm = ratio  * 449 / 12 # baseline hydro load per month
  )

# chlorophyll data, combine with loads
dat <- epcwq3 %>% 
  filter(param == "Chla") %>%
  select(date, value) %>%
  rename(chl = value) %>%
  summarise(chl = mean(chl), .by = date) |> 
  mutate(
    y = if_else(chl > 9.3, 0, 1)
  ) |> 
  inner_join(lds, by = 'date') 

TN_targets <- c(40.5, 73)
target_colors <- c("#1A99FF", "#002233")

thm <- theme_minimal() +
  theme(
    panel.grid = element_blank(), 
    axis.ticks = element_line(), 
    legend.position = 'bottom'
  )

# normalized load -----------------------------------------------------------------------------

# Fit the logistic model
mod <- glm(y ~ TN_load_norm, data = dat, family = 'binomial')

# Make predictions with confidence intervals
toprd <- seq(0, max(dat$TN_load_norm), 0.1)
plot_data <- predict(mod, newdata = tibble(TN_load_norm = toprd), type = 'response', se.fit = TRUE) %>%
  data.frame() |> 
  mutate(
    TN_load_norm = toprd,
    median = fit,
    upr = fit + 1.96 * se.fit,
    lwr = fit - 1.96 * se.fit
  )

# target predictions
# use different value for extreme scenario because outside of range
trg_data <- predict(mod, newdata = data.frame(TN_load_norm = c(TN_targets[1], max(toprd))), type = 'response', se.fit = TRUE) %>%
  data.frame() %>%
  mutate(
    grp = c('486 tons/yr', 'Largest observed'),
    TN_load_norm = c(TN_targets[1], max(toprd)),
    median = fit,
    upr = fit + 1.96 * se.fit,
    lwr = fit - 1.96 * se.fit
  )

# normalized model plot
p1 <- ggplot(plot_data, aes(x = TN_load_norm, y = median)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.3) +
  geom_line(color = "#4DB8FF", linewidth = 1) +
  geom_rect(data = trg_data, 
            aes(xmin = 0, xmax = TN_load_norm, ymin = lwr, ymax = upr, group = grp, fill = grp), 
            alpha = 0.5, color = NA
  ) +
  geom_segment(
    data = trg_data, aes(x = TN_load_norm, y = 0, xend = TN_load_norm, yend = upr, group = grp, color = grp),
    linetype = 'solid', alpha = 0.5, show.legend = FALSE
  ) + 
  scale_fill_manual(values = target_colors, name = "Loading scenario") +
  scale_color_manual(values = target_colors, name = "Loading scenario") +
  coord_cartesian(
    ylim = c(0, 1), 
    xlim = c(0, max(toprd))
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  labs(
    title = 'Probability of obtaining regulatory chlorophyll threshold',
    subtitle = "Hydrologically-normalized loading",
    x = "Tons per month",
    y = "Probability"
  ) +
  thm

# absolute load -------------------------------------------------------------------------------

# Fit the logistic model
mod <- glm(y ~ TN_load, data = dat, family = 'binomial')

# Make predictions with confidence intervals
toprd <- seq(0, 150, 0.5)
plot_data <- tibble(TN_load = toprd) %>%
  predict(mod, newdata = ., type = 'response', se.fit = TRUE) %>%
  data.frame() |> 
  mutate(
    TN_load = toprd,
    median = fit,
    upr = fit + 1.96 * se.fit,
    lwr = fit - 1.96 * se.fit
  )

# target predictions
trg_data <- predict(mod, newdata = data.frame(TN_load = TN_targets), type = 'response', se.fit = TRUE) %>%
  data.frame() %>%
  mutate(
    grp = c('486 tons/yr', 'Largest observed'),
    TN_load = TN_targets,
    median = fit,
    upr = fit + 1.96 * se.fit,
    lwr = fit - 1.96 * se.fit
  )

# Create first plot with ggplot2
p2 <- ggplot(plot_data, aes(x = TN_load, y = median)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.3) +
  geom_line(color = "#4DB8FF", linewidth = 1) +
  geom_rect(data = trg_data, 
            aes(xmin = 0, xmax = TN_load, ymin = lwr, ymax = upr, group = grp, fill = grp), 
            alpha = 0.5, color = NA
  ) +
  geom_segment(
    data = trg_data, aes(x = TN_load, y = 0, xend = TN_load, yend = upr, group = grp, color = grp),
    linetype = 'solid', alpha = 0.5, show.legend = FALSE
  ) + 
  scale_fill_manual(values = target_colors, name = "Loading scenario") +
  scale_color_manual(values = target_colors, name = "Loading scenario") +
  coord_cartesian(
    ylim = c(0, 1), 
    xlim = c(0, max(toprd))
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  labs(
    subtitle = "Absolute loading",
    x = "Tons per month",
    y = "Probability"
  ) +
  thm

# combine plots -------------------------------------------------------------------------------

p <- p1 + p2 + 
  plot_layout(ncol = 2, guides = 'collect', axis_titles = 'collect') & theme(legend.position = 'bottom')

png(here('figs/logisticmod93_revised.png'), height = 5, width = 10, units = 'in', res = 300)
print(p)
dev.off()
