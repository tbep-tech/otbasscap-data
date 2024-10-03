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
  geom_hline(yintercept = 1999.5, lwd = 1.5) + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) + 
  labs(
    title = '(a) Target (RYG management outcomes)', 
    subtitle = 'Actual data'
  )
p1b <- show_matrix(epcdatamin, bay_segment = 'OTB') + 
  geom_hline(yintercept = 1999.5, lwd = 1.5) + 
  thm + 
  labs(
    subtitle = 'Proposed targets attained'
  )
p2a <- show_wqmatrix(epcdata, param = 'chla', bay_segment = 'OTB') +
  geom_hline(yintercept = 1999.5, lwd = 1.5) + 
  thm + 
  labs(
    title = '(b) Threshold (RG regulatory outcomes)',
    subtitle = 'Actual data'
  )
p2b <- show_wqmatrix(epcdatamin, param = 'chla', bay_segment = 'OTB') + 
  geom_hline(yintercept = 1999.5, lwd = 1.5) + 
  thm + 
  labs(
    subtitle = 'Proposed targets attained'
  )

p <- p1a + p1b + p2a + p2b + plot_layout(ncol = 4)

png(here('figs/chl_hyp.png'), width = 8, height = 8, units = 'in', res = 600)
print(p)
dev.off()
