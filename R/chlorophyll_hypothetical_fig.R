library(tidyverse)
library(tbeptools)
library(patchwork)
library(here)

load(file = here::here("data/otb_subseg_sites.RData"))

nwstas <- otb_subseg_sites |> 
  filter(subsegment == 'NW') |> 
  pull(site)
cwstas <- otb_subseg_sites |> 
  filter(subsegment == 'CW') |> 
  pull(site)

# prep otb, apply min for NW/CW summer months
epcdatamin <- epcdata |> 
  mutate(
    chla = case_when(
      epchc_station %in% nwstas & mo %in% c(6:10) ~ pmin(chla, 11.3), 
      epchc_station %in% cwstas & mo %in% c(6:10) ~ pmin(chla, 13.8),
      T ~ chla
    )
  )

thm <- theme( 
  axis.text = element_blank(), 
  axis.ticks = element_blank()
)
p1a <- show_matrix(epcdata, bay_segment = 'OTB') + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) + 
  labs(
    title = '(a) Management outcomes', 
    subtitle = 'Actual data'
  )
p1b <- show_matrix(epcdatamin, bay_segment = 'OTB') + 
  thm + 
  labs(
    subtitle = 'Proposed targets attained'
  )
p2a <- show_wqmatrix(epcdata, param = 'chla', bay_segment = 'OTB') + 
  thm + 
  labs(
    title = '(b) Attainment of chlorophyll threshold',
    subtitle = 'Actual data'
  )
p2b <- show_wqmatrix(epcdatamin, param = 'chla', bay_segment = 'OTB') + 
  thm + 
  labs(
    subtitle = 'Proposed targets attained'
  )

p <- p1a + p1b + p2a + p2b + plot_layout(ncol = 4)

png(here('figs/chl_hyp.png'), width = 8, height = 8, units = 'in', res = 300)
print(p)
dev.off()
