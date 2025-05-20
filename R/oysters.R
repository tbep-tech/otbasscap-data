# setup ---------------------------------------------------------------------------------------

library(sf)
library(tidyverse)
library(here)
library(tbeptools)

trgcol <- "#1A99FF"

segs <- tbeptools::tbsegshed %>% 
  st_transform(crs = 6443)

load(file = here('data/oyse.RData'))
oyse <- oyse %>% 
  mutate(FLUCCSCODE = 6540) %>% 
  select(FLUCCSCODE, geometry = x)

# data process --------------------------------------------------------------------------------

res <- list.files(here('data'), '^sgdat') %>% 
  enframe() %>% 
  group_by(value) %>% 
  nest %>% 
  mutate(
    yr = str_extract(value, '\\d{4}')
  ) %>% 
  filter(yr >= 2014) %>% 
  mutate(
    dat = purrr::map(value, function(x){
      
      cat(x, '\t')
      
      # import file
      load(file = here(paste0('data/', x)))
      dat <- get(gsub('\\.RData', '', x))
      
      if(grepl('2022', value))
        dat <- dat %>% 
        bind_rows(oyse)
      
      dat_out <- dat |> 
        filter(FLUCCSCODE %in% c(6540)) |>
        ungroup() %>% 
        st_intersection(segs) %>% 
        summarise(
          across(geometry, st_union), 
          .by = bay_segment
        ) %>% 
        mutate(
          acre = st_area(.), 
          acre = as.numeric(units::set_units(acre, acre))
        ) %>% 
        st_set_geometry(NULL)
      
      return(dat_out)
      
    })
  ) %>% 
  ungroup() %>% 
  select(yr, dat) %>% 
  unnest(dat) %>% 
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB', 'BCB', 'MR', 'TCB'))
  )

# observed data plot --------------------------------------------------------------------------

toplo1 <- res %>% 
  arrange(bay_segment, yr) %>% 
  mutate(
    lagyr = lag(yr), 
    lagacre = lag(acre),
    .by = bay_segment
  ) %>% 
  mutate(
    chg = round(acre - lagacre, 1),
    sgn = factor(sign(chg), levels = c(-1, 0, 1), labels = c('-', '', '+')),
    chg = ifelse(is.na(sgn), NA, paste0(sgn, chg))
  ) %>% 
  filter(bay_segment == 'OTB')

p1 <- ggplot(toplo1, aes(x = yr, y = acre)) +
  geom_col(fill = '#008276', color = 'black') +
  geom_text(aes(label = chg), 
            vjust = -0.5, 
            size = 4, 
            color = 'black'
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  labs(
    x = NULL,
    y = 'Area (acres)'
  ) +
  # facet_grid(bay_segment ~ .) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

png(here("figs/oysobs.png"), width = 5, height = 4, res = 300, units = 'in')
print(p1)
dev.off()

# projected data plot -------------------------------------------------------------------------

toplo2 <- res %>% 
  mutate(
    yr = as.numeric(yr)
  ) %>% 
  filter(bay_segment == 'OTB')

mod <- lm(acre ~ yr, toplo2)
cfs <- coefficients(mod)
trgval <- 175
glyr <- (trgval - cfs[1]) / cfs[2]
toprd <- seq(2014, glyr, length.out = 100)
prds <- predict(mod, newdata = data.frame(yr = toprd)) %>% 
  tibble(acre = .) %>% 
  mutate(
    yr = toprd
  )

p2 <- ggplot(toplo2, aes(x = yr, y = acre)) +
  geom_col(fill = '#008276', color = 'black') +
  geom_line(data = prds, linetype = 'solid', color = 'darkgrey', size = 1) +
  geom_segment(
    x = 2014, xend = glyr, y = 175, yend = 175,
    color = trgcol, linewidth = 1, 
  ) +
  geom_segment(
    x = glyr, xend = glyr, y = 175, yend = 0,
    arrow = arrow(length = unit(0.75,"cm"), type = 'closed', angle = 15),
    color = trgcol, linewidth = 1, 
  ) +
  geom_text(
    x = 2014, y = trgval, label = paste(trgval, 'acres'), 
    hjust = 0, vjust = 1.2, size = 6, color = trgcol
  ) +
  geom_text(
    x = glyr, y = 0, label = round(glyr, 0), 
    hjust = 1.2, vjust = 0, size = 6, color = trgcol
  ) +
  labs(
    x = NULL,
    y = 'Area (acres)'
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

png(here("figs/oysprd.png"), width = 5, height = 4, res = 300, units = 'in')
print(p2)
dev.off()

svg(here("figs/oysprd.svg"), width = 5, height = 4, bg = 'transparent')
print(p2)
dev.off()
