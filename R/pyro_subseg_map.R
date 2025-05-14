# Required packages
library(tidyverse)
library(sf)
library(maptiles)
library(tidyterra)
library(ggspatial)
library(viridis)
library(here)

data(otbsub)
data(pyro)

# hab sum
hab_summary <- pyro %>%
  filter(!is.na(pyro)) %>%
  filter(yr <= 2021) %>%
  group_by(subsegment, yr) %>%
  summarize(max_pyro = max(pyro, na.rm = TRUE)) %>%
  group_by(subsegment) %>%
  summarize(avg_max_pyro = mean(max_pyro, na.rm = TRUE)) %>%
  mutate(avg_max_pyro = avg_max_pyro / 1e5)

# join wiht otbsub
otbsub_with_data <- otbsub %>%
  left_join(hab_summary, by = c("subseg" = "subsegment"))

bbx <- st_bbox(otbsub_with_data)
exp <- 0.01

# get basemap
basemap_tiles <- get_tiles(
  otbsub_with_data, 
  provider = 'CartoDB.PositronNoLabels',
  zoom = 12,
  crop = F
)

lbs <- st_centroid(otbsub)

# plot
p <- ggplot() +
  geom_spatraster_rgb(data = basemap_tiles) +
  geom_sf(data = otbsub_with_data, 
          aes(fill = avg_max_pyro), 
          alpha = 0.7,
          color = "darkgrey",
          size = 0.3) +
  geom_sf_text(data = lbs, aes(label = subseg), size = 5) +
  scale_fill_gradient(
    low = "#FFEBEE", 
    high = "#B71C1C",
    name = expression(10^5~"Cells/Liter"),
    labels = scales::comma
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_orienteering(
        fill = c("black", "black")
      ),
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm"),
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    style = "bar",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in")
  ) +
  coord_sf(
    expand = TRUE,
    xlim = c(bbx[1] - exp * (bbx[3] - bbx[1]),
             bbx[3] + exp * (bbx[3] - bbx[1])),
    ylim = c(bbx[2] - exp * (bbx[4] - bbx[2]),
             bbx[4] + exp * (bbx[4] - bbx[2]))
  ) +
  labs(
    title = expression(paste("Average annual max ", italic("P. bahamense"), ", 2011-2021")),
    caption = "Data source: FWC"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.ticks = element_line(),
    axis.title = element_blank(),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(1, "cm")
  )

png(here('figs/pyrosubsegmap.png'), width = 5, height = 7, units = 'in', res = 300)
print(p)
dev.off()
