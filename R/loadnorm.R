rm(list=ls(all=TRUE)) 

setwd("C:\\Users\\miles\\Documents\\GitHub\\otbasscap-data\\R")
load("../data-raw/totanndat.RData")

otbdat <- totanndat[ which( totanndat$bay_segment=="Old Tampa Bay"), ] |> as.data.frame()
head( otbdat )



# Compute hydrologically normalized TN load
hy_ref <- 449
otbdat$tn_load_norm <- otbdat$tn_load * hy_ref / otbdat$hy_load

# Plot TN loads
png( "../figs/loadnorm.png", width = 8, height = 8, units = "in", res = 600 )
par( mfrow=c(2,1), mar=c(4,5,3,5) )
plot( tn_load ~ year, data = otbdat, type = 'l',
      main = "(a) TN loads to OTB", xlab = "", ylab = "",
      ylim = c(0,1000),
      las = 1, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.3,
      col = rgb(0,0.1,0.4,0.8), lwd = 3 )

abline( v = axTicks(1), col = rgb(0,0,0,0.2) )  # grid
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )  # grid

points( tn_load ~ year, data = otbdat,  # actual load
        pch = 21, col = rgb(0,0.1,0.4,1),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

mtext( "TN load (tons)", side = 2, line = 3.5, cex = 1.3 )  # ylab
mtext( expression("TN load rate (tons /"~10^6~"m3)"),
       side = 4, line = 3.5, cex = 1.3 )  # y2 lab
axis( 4, at = seq(0,1000,200), labels = round(seq(0,2.2,length.out=6),1),     # y2
      cex.lab = 1.3, cex.axis = 1.3, las = 1)

abline( h = 486, col = rgb(0,0,0,0.5), lty = 2, lwd = 3 )  # TMDL

lines( tn_load_norm ~ year, data = otbdat,  # normalized load
       col = rgb(0.7,0.1,0.2,0.8), lwd = 3 )
points( tn_load_norm ~ year, data = otbdat, pch = 21, col = rgb(0.7,0.1,0.2,0.8),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

legend( 'bottomleft', horiz = FALSE, bty='n',
        lty = c(1,1,2), lwd = 3, cex = 1.2,
        col = c( rgb(0,0.1,0.4,0.8), rgb(0.7,0.1,0.2,0.8), rgb(0,0,0,0.5) ),
        legend = c("Actual load","Normalized load",
                   expression("TMDL (486 tons/yr) = Delivery ratio (1.08 tons /"*10^6~"m3)") )
         )


# Plot hydro loads
plot( hy_load ~ year, data = otbdat, type = 'l',
      main = "(b) Hydrologic loads to OTB", xlab = "", ylab = "",
      las = 1, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.3,
      col = rgb(0,0.4,0.7,0.8), lwd = 3 )

abline( v = axTicks(1), col = rgb(0,0,0,0.2) )  # grid
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )  # grid

abline( h = hy_ref, col = rgb(0,0,0,0.5), lty = 2, lwd = 3 )  # reference hydrologic load
text( x = 2012, y = hy_ref, pos = 1, col = rgb(0,0,0,0.7),
      labels = "Reference hydrologic load (449 million m3)" )

points( hy_load ~ year, data = otbdat, pch = 21, col = rgb(0,0.4,0.7,1),
        cex = 1.2, lwd = 2, bg = rgb(1,1,1,1) )

mtext( "million m3", side = 2, line = 3.5, cex = 1.3 )  # ylab



dev.off()

# recreate with ggplot ------------------------------------------------------------------------

library(tidyverse)
library(here)

load(url('https://github.com/tbep-tech/load-estimates/raw/refs/heads/main/data/totanndat.RData'))

hy_ref <- 449
tn_ref <- 468
ptsz <- 2.5
  
otbdat <- totanndat |> 
  filter(bay_segment == 'Old Tampa Bay') |> 
  mutate(
    tn_load_norm = tn_load * hy_ref / hy_load
  )

toplo1 <- otbdat |> 
  select(year, tn_load, tn_load_norm) |> 
  pivot_longer(cols = c(tn_load, tn_load_norm), names_to = 'load_type', values_to = 'load') |> 
  mutate(
    load_type = ifelse(load_type == 'tn_load', 'Actual', 'Normalized')
  )

# ggplot with transparent gray background from 1992 to 1994
p1 <- ggplot() + 
  geom_rect(data = data.frame(x = 0, y = 0), aes(xmin = 1991.5, xmax = 1994.5, ymin = -Inf, ymax = Inf, fill = 'Reference period'), color = NA, alpha = 0.6) +
  geom_hline(aes(yintercept = tn_ref, linetype = paste0('Threshold (', tn_ref, ' tons / yr)')), color = 'black') +
  geom_line(data = toplo1, aes(x = year, y = load, color = load_type), linewidth = 1) + 
  geom_point(data = toplo1, aes(x = year, y = load, color = load_type), pch = 21, fill = 'white', size = ptsz, stroke = 2) + 
  scale_color_manual(values = c(rgb(0,0.1,0.4,1), rgb(0.7,0.1,0.2,1))) +
  scale_fill_manual(values = 'gray') +
  scale_linetype_manual(values = 'dashed') +
  theme_minimal() + 
  theme(
    legend.position = 'top',
    panel.grid.minor = element_blank()
  ) + 
  labs(
    x = NULL,
    linetype = NULL, 
    color = NULL, 
    y = expression(paste('tons'~yr^{-1})),
    title = 'Nitrogen loads',
    fill = NULL
  )

p2 <- ggplot() + 
  geom_rect(data = data.frame(x = 0, y = 0), aes(xmin = 1991.5, xmax = 1994.5, ymin = -Inf, ymax = Inf, fill = 'Reference period'), color = NA, alpha = 0.6) +
  geom_hline(aes(yintercept = hy_ref, linetype = paste0('Reference hydrologic load (', hy_ref, ' million m3)')), color = 'black') +
  geom_line(data = otbdat, aes(x = year, y = hy_load), linewidth = 1, color = rgb(0,0.4,0.7,1)) + 
  geom_point(data = otbdat, aes(x = year, y = hy_load), pch = 21, fill = 'white', color = rgb(0,0.4,0.7,1), size = ptsz, stroke = 2) + 
  scale_fill_manual(values = 'gray') +
  scale_linetype_manual(values = 'dashed') +
  theme_minimal() + 
  theme(
    legend.position = 'top',
    panel.grid.minor = element_blank()
  ) + 
  labs(
    x = NULL,
    linetype = NULL, 
    title = 'Hydrologic loads',
    y = expression(paste(10^6~m^3~yr^{-1})),
    fill = NULL
  )

png(here('figs/loadnorm1withref.png'), width = 7, height = 3, units = 'in', res = 300)
print(p1)
dev.off()

svg(here('figs/loadnorm1withref.svg'), width = 7, height = 3, bg = 'transparent')
print(p1)
dev.off()

png(here('figs/loadnorm2withref.png'), width = 7, height = 3, units = 'in', res = 300)
print(p2)
dev.off()

svg(here('figs/loadnorm2withref.svg'), width = 7, height = 3, bg = 'transparent')
print(p2)
dev.off()


