rm(list=ls(all=TRUE))

# Load signals
load("../data/ssa_TN_otb.RData")
load("../data/ssa_TP_otb.RData")
load("../data/ssa_sal_otb.RData")
load("../data/ssa_temp_otb.RData")
load("../data/ssa_Qtarpon_otb.RData")

# Generate figure
obj <- ls()[c(4,5,3,2,1)]
objlab <- c( "TN load", "TP load", "Temperature", "Salinity","Lake Tarpon canal discharge")
units <- c( "mg/L", "mg/L", "C", "ppt",
            expression(10^9~"ft3/mo") )

png( "../figs/ssa_driver_plots.png", height = 7.5, width = 7, units = 'in', res = 500 )
  par( mfrow=c(5,1), mar=c(1,5,2,1), oma=c(2,0,0,0) )
  for( i in 1:length( obj ) ){
    
    this <- eval( parse(text=obj[i]) )
    
    plot( ts ~ month, dat = this$dat,
          type = 'l', las = 1,
          xlab = "", ylab = "", yaxt = 'n', xaxt = 'n',
          lwd = 2, col = rgb(0,0,0,0.6) )
    mtext( paste0( objlab[i],"  (signal strength: ",
                   round(this$sigstrength,2), "%)" ),
           side = 3, line = 0,
           font = 1, cex = 1.1, adj = 0 )
    date.seq <- seq.Date( min(this$dat$month), max(this$dat$month)+365, 'year' )
    if( i==length(obj) ){
      axis( 1, at = date.seq, labels = year(date.seq), cex.axis = 1.3  )
    } else {
      axis( 1, at = date.seq, labels = rep("",length(date.seq))  )
    }
    if( obj[i] == "ssa_Qtarpon_otb" ){
      axis( 2, at = axTicks(2), labels = axTicks(2)/1e9, las = 1, cex.axis = 1.2 ) 
    } else {
      axis( 2, las = 1, cex.axis = 1.2 )
    }
    mtext( units[i],
           side = 2, line = 3, cex = 1 )
    lines( signal ~ month, dat = this$dat,
           lwd = 3, col = rgb(1,0.2,0.1,0.8) )
    abline( v = date.seq, col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    
    if( i==1 ){
      legend( "topright", horiz = TRUE, bty = 'n',
              xpd = TRUE, inset = c(0,-0.35),
              legend = c("data","signal"), lwd = c(2,4), cex = 1.1,
              col = c( rgb(0,0,0,0.6), rgb(1,0.2,0.1,0.8) ) )
    }
    
  }  # // end i loop
  
dev.off()
