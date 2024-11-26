library(haven)
library(dplyr)
library(sf)
library(leaflet)
library(tibble)
library(tidyr)
library(leaflet)
library(mapedit)

load(file = here::here('data/otbsub.RData'))

pth <- 'T:/03_BOARDS_COMMITTEES/05_TBNMC/TB_LOADS/2022_RA_Deliverables/2017-2021Annual&MonthlyLoadDatasets/NPS1721/Monthly'

dat <- read_sas(file.path(pth, 'allnuts1721monthbasin20221027.sas7bdat')) |> 
  filter(BAY_SEG == 1)

shp <- st_read('T:/05_GIS/BOUNDARIES/TBEP_dBasins_Correct_Projection.shp') |> 
  filter(BAYSEGNAME == 'Old Tampa Bay') |> 
  st_transform(crs = 4326)

gagebas <- shp |> 
  mutate(
    subbasin = case_when(
      NEWGAGE %in% 'LTARPON' ~ 'NW', 
      NEWGAGE %in% '02307359' ~ 'NW', 
      NEWGAGE %in% c('02307000', '02306647') ~ 'NE', 
      T ~ NA_character_   
    )
  ) |> 
  filter(!is.na(subbasin))

# interactively select polygons on mapiew

# create leaflet polygon map with ungage
m <- leaflet() |> 
  addProviderTiles('CartoDB.Positron') |> 
  addPolygons(data = shp[shp$NEWGAGE == '206-1', ], weight =  0.5) |> 
  addPolygons(data = otbsub) |> 
  addDrawToolbar(
    targetGroup = 'draw',
    editOptions = editToolbarOptions()
  )

NW <- selectFeatures(ungage, map = m)
NE <- selectFeatures(ungage, map = m)
CW <- selectFeatures(ungage, map = m)
CE <- selectFeatures(ungage, map = m)
SW <- selectFeatures(ungage, map = m)
SE <- selectFeatures(ungage, map = m)

ungagebas <- list(
  NW = NW,  
  NE = NE,
  CW = CW,
  CE = CE,
  SW = SW,
  SE = SE
  ) |> 
  enframe(name = 'subbasin', value = 'data') |> 
  unnest('data') |> 
  st_sf()

allsubbas <- bind_rows(ungagebas, gagebas) 

save(allsubbas, file = here::here('data/allsubbas.RData'))

tmp <- allsubbas |> 
  group_by(NEWGAGE, subbasin) |> 
  summarise() |> 
  ungroup()
tmp$acres <- as.numeric(st_area(tmp) / 4047)
tmp <- tmp |> 
  mutate(
    prop = case_when(
      NEWGAGE != '206-1' ~ 1,
      TRUE ~ acres / sum(acres)
    )
  ) |> 
  st_set_geometry(NULL)

# tmp2 <- dat |> 
#   select(basin, source, year, month, tnload) |> 
#   left_join(tmp, by = c('basin' = 'NEWGAGE'), relationship = 'many-to-many') |> 
#   mutate(tnload = tnload * prop)
# # verify proportional sums are same as dat sums for tn load
