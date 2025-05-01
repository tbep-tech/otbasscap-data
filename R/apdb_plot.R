# setup ---------------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(patchwork)
library(here)

# data prep -----------------------------------------------------------------------------------

datraw <- read_excel(here('data-raw/DUMP_apdb.xlsx'))

dat <- datraw |> 
  select(yr = CompletionDate, ProjectCompleted, SegmentName, TN_lbs_yr) |> 
  filter(ProjectCompleted == 1) |> 
  filter(yr != 0) |> 
  filter(yr >= 2002 & yr <= 2026) |> 
  mutate(
    raperiod = case_when(
      yr <= 2026 & yr >= 2022 ~ '2022-2026',
      yr <= 2021 & yr >= 2017 ~ '2017-2021',
      yr <= 2016 & yr >= 2012 ~ '2012-2016',
      yr <= 2011 & yr >= 2007 ~ '2007-2011',
      yr <= 2006 & yr >= 2002 ~ '2002-2006'
    ), 
    raperiod = factor(raperiod, levels = rev(c('2022-2026', '2017-2021', '2012-2016', '2007-2011', '2002-2006')))
  ) |> 
  filter(SegmentName %in% c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay', 'Boca Ciega Bay', 'Terra Ceia Bay', 'Manatee River')) |> 
  mutate(
    SegmentName = factor(SegmentName, levels = c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay', 'Boca Ciega Bay', 'Terra Ceia Bay', 'Manatee River'))
  )

toplo <- dat |> 
  summarise(
    TN_lbs_yr = sum(TN_lbs_yr, na.rm = T),
    cnt = n(),
    .by = c('raperiod', 'SegmentName')
  ) |> 
  mutate(
    TNtons_period = 5 * TN_lbs_yr / 2000,
  )

# all bay segments ----------------------------------------------------------------------------

cols <- c("#00577C", "#ABB8CE", "#077686", "#DA8C53", "#4E795E", "#F1BCB5", "#FDF486")

# ggplot of cnt by ra period and bay segment
p1 <- ggplot(toplo, aes(x = raperiod, y = cnt, fill = SegmentName)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = cols) +
  labs(
    title = 'Reported Projects by Bay Segment and Reporting Period',
    subtitle = 'Number of Projects',
    x = NULL,
    y = 'Count'
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank()
  )

# ggplot of TNtons_period by ra period and bay segment
p2 <- ggplot(toplo, aes(x = raperiod, y = TNtons_period, fill = SegmentName)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = cols) +
  labs(
    subtitle = 'TN Tons Removed',
    x = 'RA Period',
    y = 'Tons', 
    caption = 'Note: Includes projects with a completion year and reported TN load reduction'
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank()
  )

p <- p1 + p2 + guide_area() + plot_layout(ncol = 1, guides = 'collect', heights = c(1, 1, 0.4))

png(here('figs', 'apdb_tampa_bay.png'), width = 8, height = 6, units = 'in', res = 300)
print(p)
dev.off()

# old tampa bay only --------------------------------------------------------------------------

toplo2 <- toplo |> 
  filter(SegmentName == 'Old Tampa Bay')

# ggplot of cnt by ra period and bay segment
p1 <- ggplot(toplo2, aes(x = raperiod, y = cnt)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(
    title = 'Reported Projects and Reporting Period, Old Tampa Bay only',
    subtitle = 'Number of Projects',
    x = NULL,
    y = 'Count'
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank()
  )

# ggplot of TNtons_period by ra period and bay segment
p2 <- ggplot(toplo2, aes(x = raperiod, y = TNtons_period)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(
    subtitle = 'TN Tons Removed',
    x = 'RA Period',
    y = 'Tons', 
    caption = 'Note: Includes projects with a completion year and reported TN load reduction'
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank()
  )

p <- p1 + p2 + plot_layout(ncol = 1)

png(here('figs', 'apdb_old_tampa_bay.png'), width = 8, height = 5, units = 'in', res = 300)
print(p)
dev.off()
