# Required packages
library(tidyverse)
library(sf)
library(maptiles)
library(tidyterra)
library(ggspatial)
library(viridis)
library(here)
library(tbeptools)

data(otbsub)

chldat <- epcdata |> 
  anlz_avedatsite() |> 
  anlz_attainsite(thr = 'chla') |> 
  filter(bay_segment == 'OTB') |> 
  filter(yr >= 2011 & yr <= 2021) |> 
  summarise(
    val = mean(val, na.rm = T),
    .by = epchc_station
  ) |> 
  mutate(
    met = ifelse(val < 9.3, 'met', 'not met')
  ) |> 
  left_join(stations, by = 'epchc_station') |> 
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326)
  
bbx <- st_bbox(otbsub)
exp <- 0.01

# get basemap
basemap_tiles <- get_tiles(
  otbsub, 
  provider = 'CartoDB.PositronNoLabels',
  zoom = 12,
  crop = F
)

# plot
p <- ggplot() +
  geom_spatraster_rgb(data = basemap_tiles) +
  geom_sf(data = otbsub, alpha = 0.7, fill = NA, color = "black") +
  geom_sf(data = chldat, 
          aes(color = met, size = val)
  ) +
  scale_colour_manual(values = c('#2DC938', '#CC3231')) +
  guides(size = 'none', color = guide_legend(override.aes = list(size = 5))) +
  ggrepel::geom_text_repel(data = chldat, aes(label = round(val, 1), geometry = geometry), stat = "sf_coordinates", size = 3, inherit.aes = F) +
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
    title = "Mean annual avg chl-a (ug/L), 2011-2021",
    caption = "Data source: EPCHC", 
    color = "Threshold"
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

png(here('figs/chlasubsegmap.png'), width = 4.5, height = 6, units = 'in', res = 300)
print(p)
dev.off()
