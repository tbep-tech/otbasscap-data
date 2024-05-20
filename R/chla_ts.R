rm(list=ls(all=TRUE)) 

load(file = here::here('data-clean/epcwq_clean.RData'))

subsegs <- c("NW","NE","CW","CE","SW","SE")

par(mfrow=c(3,2), mar=c(4,4,3,1))

# Plot data and annual averages over time at each OTB subsegment
for( i in 1:length(subsegs) ){
  
  this.dat <- epcwq3[ which( epcwq3$subseg==subsegs[i] &
                               epcwq3$param == "Chla" ),]
  
  # plot monthly data
  plot( value ~ date, type = 'l',
        data = this.dat, ylim = c(0,30), las = 1,
        main = paste0("Chla at OTB ",subsegs[i]," segment"),
        xlab = "", ylab = "Chla (ug/L)",
        pch = 16, cex = 0.7, col = rgb(0,0,0.7,0.2) )
  # lines( value ~ date,
  #        data = this.dat, col = rgb(0,0,0.3,0.2) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  abline( v = seq.Date( as.Date("2000-01-01"), as.Date("2024-01-01"), 'year'),
          col = rgb(0,0,0,0.1) )
  
  # draw target
  abline( h = 9.3, lty = 2, lwd = 2, col = rgb(0.3,0.1,0,0.8) )
  
  # plot annual means
  this.dat$year <- this.dat$date |> year()
  means <- this.dat |> dplyr::summarise( value = mean(value),
                                         .by = year )
  means$date <- paste0(means$year,"-07-01") |> as.Date()
  sd <- this.dat |> dplyr::summarise( sd = sd(value),
                                      .by = year )
  n <- this.dat |> dplyr::summarise( n = length(value),
                                     .by = year )
  lines( value ~ date, data = means,
         lwd = 2, col = rgb(0.8,0.1,0,0.5) )
  segments( x0 = means$date,
            y0 = c( means$value - sd$sd/sqrt(n$n) ),
            y1 = c( means$value + sd$sd/sqrt(n$n) ),
            lwd = 2, col = rgb(0.8,0.1,0,0.8)
  )
  points( value ~ date, data = means,
          pch = 21, cex = 1.3, lwd = 1.5,
          col = rgb(0.8,0.1,0,1), bg = rgb(1,1,1,1) )
  
}


# Plot data by month-of-year with monthly averages at each OTB subsegment
for( i in 1:length(subsegs) ){
  
  this.dat <- epcwq3[ which( epcwq3$subseg==subsegs[i] &
                               epcwq3$param == "Chla" ),]
  # this.dat$year <- this.dat$date |> year()
  
  # initialize plot
  plot( x=0, y=0, xaxt = 'n',
        xlim = c(1,12), ylim = c(0,50), las = 1,
        main = paste0("Chla at OTB ",subsegs[i]," segment"),
        xlab = "", ylab = "Chla (ug/L)",
        pch = 16, cex = 0.7, col = rgb(0,0,0,0) )
  # lines( value ~ date,
  #        data = this.dat, col = rgb(0,0,0.3,0.2) )
  axis( 1, at = 1:12, labels = month(1:12,abbr=TRUE,label=TRUE) )
  abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
  abline( v = 1:12, col = rgb(0,0,0,0.1) )
  
  for( j in unique(year(this.dat$date)) ){
    # monthly data
    this.dat.year <- this.dat[ which( year(this.dat$date) == j ), ]
    this.dat.year$month <- this.dat.year$date |> month()
    this.n <- nrow(this.dat.year)
    points( x = this.dat.year$month + rnorm(this.n,0,0.05),
            y = this.dat.year$value,
            pch = 16, cex = 0.7,
           col = rgb( 0.1, 0, ((j-2000+2)/30), ((j-2000+2)/50) ) )
    
  }  # // end j loop

  # draw target
  abline( h = 9.3, lty = 2, lwd = 2, col = rgb(0.3,0.1,0,0.8) )
    
  # plot monthly means
  this.dat$month <- this.dat$date |> month()
  means <- this.dat |> dplyr::summarise( value = mean(value),
                                         .by = month )
  sd <- this.dat |> dplyr::summarise( sd = sd(value),
                                      .by = month )
  n <- this.dat |> dplyr::summarise( n = length(value),
                                     .by = month )
  lines( value ~ month, data = means,
         lwd = 2, col = rgb(0.8,0.1,0,0.5) )
  segments( x0 = means$month,
            y0 = c( means$value - sd$sd/sqrt(n$n) ),
            y1 = c( means$value + sd$sd/sqrt(n$n) ),
            lwd = 2, col = rgb(0.8,0.1,0,0.8)
            )
  points( value ~ month, data = means,
          pch = 21, cex = 1.3, lwd = 1.5,
          col = rgb(0.8,0.1,0,1), bg = rgb(1,1,1,1) )

} # // end i loop



