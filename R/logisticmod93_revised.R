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

# prep chl
chl <- epcwq3 %>% 
  filter(param == "Chla") %>%
  select(date, value) %>%
  rename(chl = value) %>%
  summarise(chl = mean(chl), .by = date) |> 
  mutate(
    y = if_else(chl > 9.3, 0, 1)
  )

# hydro loads
hyd <- hydro %>%
  filter(bay_segment == "Old Tampa Bay", year >= 2000) %>%
  mutate(date = as.Date(paste0(year, "-", month, "-01"))) %>%
  select(year, date, hy_load_106_m3_mo) %>%
  rename(hydro_load = hy_load_106_m3_mo) |> 
  mutate(
    hydroann_load = sum(hydro_load), 
    .by = year
  )

# prep loading data
lds <- loads %>%
  filter(param == "TN load", year(date) >= 2000) %>%
  select(date, value) %>%
  rename(TN_load = value) |> 
  inner_join(hyd, by = 'date') |> 
  mutate(
    ratio = TN_load / hydro_load,
    TN_load_norm = ratio  * 449 / 12 # baseline hydro load per month
  )

# combine data
dat <- inner_join(chl, lds, by = 'date') 

# normalized load -----------------------------------------------------------------------------

# Fit the logistic model
mod <- glm(y ~ ratio, data = dat, family = 'binomial')

# Generate predictions
toprd <- tibble(ratio = seq(0, max(dat$ratio), 0.01))

# Make predictions with confidence intervals
preds <- predict(mod, newdata = toprd, type = 'response', se.fit = TRUE)
plot_data <- toprd %>%
  mutate(
    median = preds$fit,
    upr = preds$fit + 1.96 * preds$se.fit,
    lwr = preds$fit - 1.96 * preds$se.fit
  )

# Define delivery ratio targets
deliv_targets <- c(1.08, 1.71)
target_labels <- paste0(deliv_targets, " tons/Mm3")

# Create second plot with ggplot2
p1 <- ggplot(plot_data, aes(x = ratio, y = median)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.3) +
  geom_line(color = "#4DB8FF", size = 1) +
  geom_vline(xintercept = deliv_targets, linetype = "dashed", color = target_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_continuous(
    limits = c(0, 1.7),
    breaks = seq(0, 1.7, 0.2)
  ) +
  labs(
    title = "(b) Chl-a threshold attainment and TN delivery ratio",
    x = "TN delivery ratio (tons/Mm3 per year)",
    y = "Probability"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0, size = 12),
    axis.title = element_text(size = 10)
  )

# absolute load -------------------------------------------------------------------------------

# Fit the logistic model
mod <- glm(y ~ TN_load, data = dat, family = 'binomial')

# Generate predictions
toprd <- tibble(TN_load = seq(0, max(dat$TN_load), 0.5))

# Make predictions with confidence intervals
preds <- predict(mod, newdata = toprd, type = 'response', se.fit = TRUE)
plot_data <- toprd %>%
  mutate(
    median = preds$fit,
    upr = preds$fit + 1.96 * preds$se.fit,
    lwr = preds$fit - 1.96 * preds$se.fit
  )

# Define TN targets
TN_targets <- c(40.5, 73)
target_labels <- paste0(TN_targets, " tons/month")
target_colors <- c("#1A99FF", "#002233")

trgpreds <- predict(mod, newdata = data.frame(TN_load = TN_targets), type = 'response', se.fit = TRUE)
trg_data <- data.frame(
  grp = c('486 tons/yr', 'Largest observed'),
  TN_load = TN_targets,
  median = trgpreds$fit,
  upr = trgpreds$fit + 1.96 * trgpreds$se.fit,
  lwr = trgpreds$fit - 1.96 * trgpreds$se.fit
)

# Create first plot with ggplot2
p2 <- ggplot(plot_data, aes(x = TN_load, y = median)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.3) +
  geom_line(color = "#4DB8FF", size = 1) +
  geom_rect(data = trg_data, 
            aes(xmin = 0, xmax = TN_load, ymin = lwr, ymax = upr, group = grp, fill = grp), 
            alpha = 0.5, color = NA
  ) +
  geom_segment(
    data = trg_data, aes(x = TN_load, y = 0, xend = TN_load, yend = upr, group = grp, color = grp),
    linetype = 'dashed', alpha = 0.5, show.legend = FALSE
  ) + 
  scale_fill_manual(values = target_colors, name = "Loading scenario") +
  scale_color_manual(values = target_colors, name = "Loading scenario") +
  coord_cartesian(
    ylim = c(0, 1), 
    xlim = c(0, 150)
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  labs(
    subtitle = "Absolute loading",
    x = "TN load (tons/month)",
    y = "Probability"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), 
    axis.ticks = element_line()
  )


# Combine both plots using patchwork
combined_plot <- p1 / p2
