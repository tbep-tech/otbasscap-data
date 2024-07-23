rm(list=ls(all=TRUE))

#  Load libraries
if(!require(lubridate)) { install.packages('lubridate') }; library(lubridate)
if(!require(dplyr)) { install.packages('dplyr') }; library(dplyr)
if(!require(rEDM)) { install.packages('rEDM') }; library(rEDM)


# Input data --------------------------------------------------------------

# Load embedded signal objects
sigpaths <- dir("../data")[ grep("embed_",dir("../data")) ]
sigpaths <- paste0("../data/",sigpaths)
for( i in 1:length(sigpaths) ){
  load( sigpaths[i] )
}  # // end loading loop

# Assemble signals for CCM analyses
for( i in 1:length(sigpaths) ){
  objname <- sub("../data/","",sigpaths[i])
  objname <- sub(".RData","",objname)
  varname <- sub("embed_","",objname)
  obj <- eval(parse(text=objname))
  if(i==1){
    signals <- data.frame( month = obj$t,
                           sig = obj$x )
    colnames(signals)[i+1] <- varname
  } else {
    signals <- full_join( signals,
                          data.frame( month = obj$t, sig = obj$x ),
                          'month' )
    colnames(signals)[i+1] <- varname
    signals <- signals[ order(signals$month), ]
  }
}  # end signal assembly loop
colnames( signals )[1] <- "date" 


# Assemble embedding parameters for response variable (chl)
embnames <- ls()[ grep("embed_chl",ls()) ]  # get obj names
embed_tbl <- data.frame( var = NA, d = NA, tw = NA, m = NA )  # init dataframe
for( i in 1:length(embnames) ){
  obj <- eval(parse(text=embnames[i]))
  embed_tbl[i,] <- rep(NA,ncol(embed_tbl))
  embed_tbl$var[i] <- sub("embed_","",embnames[i])
  embed_tbl$d[i] <- obj$d
  embed_tbl$tw[i] <- obj$tw
  embed_tbl$m[i] <- obj$m
}


# CCM parameters ----------------------------------------------------------

delay.seq <- seq(-7,7,1)   # delay sequence for Extended CCM (absolute value of 'from' should be >= 'to' value )
min.date  <- as.Date( '2012-01-01')
max.date  <- as.Date( '2021-12-31' )
lib.step  <- 10     # library step size
ccm.frac  <- 0.10  # xmap skill will be estimated from the largest 'ccm.frac' of library sizes
nsurr     <- 500   # Number of surrogates for significance testing
write.png <- TRUE  # write plots to file (PNG)?
gpath     <- "../figs/"  # directory path for exporting images


# CCM functions -----------------------------------------------------------

#  Define function to perform bidirectional extended CCM with significance testing.
##  Many arguments default to parameters specified at the top of this script. Arguments:
##     DATA        dataframe containing "date" as first column and input signals in subsequent columns
##     LIBRARY     character string containing column name of library variable (LIBRARY xmap TARGET)
##     TARGET      character string containing column name of target variable (LIBRARY xmap TARGET)
##     EMB         embedding parameters; list containing d, tw, and m values
##     DELAY.SEQ   numeric vector containing delays for extended CCM (symmetric about zero)
##     LIB.STEP    library step size (integer)
##     SAMPLES     Number of replicate CCM tests to run (integer); defaults to 100
##     DATE.DIFF   data time step (in days); the code checks that there are no missing values. Use NULL to skip this check 
##     SEED        random seed
delay.CCM <- function( DATA=signals, LIBRARY, TARGET, EMB,
                       DELAY.SEQ=seq(delay.seq[1],-delay.seq[1],1),
                       LIB.STEP=lib.step, SAMPLES=100,
                       DATES=c(min.date,max.date), DATE.DIFF=NULL, SEED=78152937 ){
  
  start.time <- Sys.time()  # record starting time of execution
  
  # Generate variable names, replacing periods with spaces
  var.response <- gsub( '\\.', ' ', LIBRARY )
  var.driver   <- gsub( '\\.', ' ', TARGET )
  
  # Assemble data
  # Isolate date, library, and signal columns
  this.df <- DATA[,c("date",LIBRARY,TARGET)] |> na.omit()
  # Subset data by date
  if( !is.null(DATES) ){
    this.df <- this.df[ which( this.df$date >= DATES[1] & this.df$date <= DATES[2] ) ,]
  }  # // end if(DATES)
  
  # Check for missing values
  if( !is.null(DATE.DIFF) ){
    date.diff <- this.df$date |> diff() |> unique()
    if( date.diff != DATE.DIFF ){ stop('Unexpected date diffs. Check for missing data in library and target.') }
  }  # // end if(DATE.DIFF)
  # Record date range
  date.range <- range( this.df$date )
  
  # Run extended CCM at each specified delay
  # Initialize objects to hold results at each delay
  ccm.xy <- vector( mode='list', length=length(DELAY.SEQ) )  # empty list for lib:target CCM output
  names( ccm.xy ) <- paste0('delay',DELAY.SEQ)  # name list elements using delay values
  ccm.yx <- vector( mode='list', length=length(DELAY.SEQ) )  # empty list for target:lib CCM output
  names( ccm.yx ) <- paste0('delay',DELAY.SEQ)  # name list elements using delay values
  corr <- c()  # empty vector for Pearson correlations
  # Get embedding parameters from library vector
  d <- EMB$d
  tw <- EMB$tw
  m <- EMB$m
  # Run CCM at each delay
  for( i in DELAY.SEQ ){
    
    # Initialize this.data dataframe
    this.data <- data.frame()
    
    # Adjust library and target vectors for delay
    if( i<0 ){  # Negative delay: Response (library) predicts past values of driver (target)
      this.lib    <- this.df[ (-i+1):nrow(this.df) , 2 ] # for delay=-1, lib takes 2nd row (-(-1)+1=2) thru final row..
      this.target <- this.df[ 1:(nrow(this.df)+i)  , 3 ] # ..and target takes 1st row thru penultimate row (nrow+(-1)).
      this.data   <- data.frame( t=1:length(this.lib), lib=this.lib, target=this.target )
    } else if ( i>0 ){  # Positive delay: Response (library) predicts future values of driver (target)
      this.lib    <- this.df[ 1:(nrow(this.df)-i) , 2 ] # for delay=1, lib takes first row through penultimate row..
      this.target <- this.df[ (i+1):nrow(this.df) , 3 ] # ..and target takes 2nd row (1+1=2) through final row.
      this.data   <- data.frame( t=1:length(this.lib), lib = this.lib, target = this.target )
    } else {  # No delay: Response predicts concurrent values of driver
      this.data <- data.frame( t=1:nrow(this.df), lib = this.df[,2], target = this.df[,3] )
    }  # // end if(delay)
    # Compute Pearson correlation
    corr[ which(DELAY.SEQ==i) ] <- cor( this.data$lib, this.data$target ) 
    # Scale the library and target vectors to zero mean and unit variance
    this.data[,2:3] <- scale( this.data[,2:3] )
    # Generate library size sequence
    lib.max <- nrow(this.data) - (m-1)*d  # with replacement=FALSE in CCM() call below, max library size is constrained by this formula
    lib.seq <- seq( lib.max%%LIB.STEP, lib.max, LIB.STEP )
    if( lib.seq[1]<LIB.STEP ){ lib.seq <- lib.seq[-1] }
    # lib.seq <- seq(LIB.STEP,lib.max,LIB.STEP)
    # lib.seq.naive <- seq(LIB.STEP,nrow(this.data),LIB.STEP)
    # seq.diff <- nrow(this.data) - max(lib.seq.naive) # difference between sequence max and data length
    # lib.seq <- lib.seq.naive + seq.diff  # shift sequence so that sequence max equals data length
    # Run CCM on this.data
    this.ccm <- CCM( dataFrame = this.data, E = m, tau = -d, Tp = 0, exclusionRadius = tw,
                     columns = 'lib', target = 'target',
                     libSizes = lib.seq, sample = SAMPLES,
                     seed = (SEED+i), includeData = TRUE, showPlot = FALSE )
    # this.ccm <- ccm( block = this.data, 
    #                  E = m, tau = -d, exclusion_radius = tw,
    #                  num_neighbors = 'E+1', lib_sizes = c(LIB.STEP,lib.max,LIB.STEP), #lib.seq, #c(LIB.STEP,nrow(this.data),LIB.STEP),
    #                  num_samples = 100, replace = FALSE, stats_only = FALSE,
    #                  lib_column = 1, target_column = 2, RNGseed = (SEED+i) )
    # Format CCM output: ccm.xy contains lib:target results; ccm.yx contains target:lib results.
    # lib:target 
    ccm.xy[[ which( DELAY.SEQ ==  i ) ]] <- list( LibMeans = this.ccm$LibMeans,
                                                  CCM1_PredictStat = this.ccm$CCM1_PredictStat,
                                                  CCM1_Predictions = this.ccm$CCM1_Predictions )
    # Remove target:lib column from LibMeans dataframe
    ccm.xy[[ which( DELAY.SEQ == i ) ]]$LibMeans <- ccm.xy[[ which( DELAY.SEQ == i ) ]]$LibMeans[,-3]
    # target:lib
    ccm.yx[[ which( DELAY.SEQ == -i ) ]] <- list( LibMeans = this.ccm$LibMeans,
                                                  CCM2_PredictStat = this.ccm$CCM2_PredictStat,
                                                  CCM2_Predictions = this.ccm$CCM2_Predictions )
    # Remove lib:target column from LibMeans dataframe
    ccm.yx[[ which( DELAY.SEQ == -i ) ]]$LibMeans <- ccm.yx[[ which( DELAY.SEQ == -i ) ]]$LibMeans[,-2]
    
  }  # // end i loop(DELAY.SEQ)
  
  print( difftime(Sys.time(), start.time) )  # Print execution time to console
  
  # Output
  return( list( data = this.data, ccm.xy = ccm.xy, ccm.yx = ccm.yx,
                var.response = var.response, var.driver = var.driver, dates = date.range,
                lib.seq = lib.seq, corr = corr, d = d, tw = tw, m = m ) )
  # return( list( ccm.out = ccm.out, rho = rho, rho.sd = rho.sd, date.range = date.range ) )
  
}  # // end delay.CCM()



#  Define function to plot output from delay.CCM() and compute rho estimates
##  Most arguments default to parameters specified at the top of this script. Arguments:
##     CCM.OUTPUT  list object containing delay.CCM() output
##     MEAN.FRAC   the fraction of libraries from which to compute mean and sd rho values (xmap skill)
##     DELAY.SEQ   numeric vector containing delays for extended CCM
##     CORR        logical indicating whether to plot Pearson correlation on rho vs. library size plots
##     WRITE.PNG   logical indicating whether to write plots to file (PNG)
plot.CCM <- function( CCM.OUTPUT, MEAN.FRAC=ccm.frac, DELAY.SEQ=delay.seq,
                      PLOT.CORR=TRUE, WRITE.PNG=write.png ){
  
  # Load CCM output
  ccm.xy       <- CCM.OUTPUT$ccm.xy
  ccm.yx       <- CCM.OUTPUT$ccm.yx
  var.response <- CCM.OUTPUT$var.response
  var.driver   <- CCM.OUTPUT$var.driver
  corr         <- CCM.OUTPUT$corr
  dates        <- CCM.OUTPUT$dates
  
  # Initialize vectors to hold rho means and sd's
  rho.xy  <- rep(NA,length(DELAY.SEQ))  # mean rho values for lib:target
  rho.sd.xy <- rep(NA,length(DELAY.SEQ))  # rho standard deviation values for lib:target
  rho.yx  <- rep(NA,length(DELAY.SEQ))  # mean rho values for lib:target
  rho.sd.yx <- rep(NA,length(DELAY.SEQ))  # rho standard deviation values for lib:target
  
  # Plot rho vs. library size at each non-positive delay
  if(write.png){ 
    png( paste0(gpath,'ccm__',var.response,'_X_',var.driver,'__',dates[1],'--',dates[2],'__rho-vs-lib.png'),
         height=11, width=9, res=600, units='in' ) }  # // end if(write.png)
  par(mfrow=c(4,2))
  for( i in DELAY.SEQ ){
    # Get results for the delay i
    this.ccm.xy <- ccm.xy[[ which( DELAY.SEQ == i ) ]]
    this.ccm.yx <- ccm.yx[[ which( DELAY.SEQ == i ) ]]
    lib.seq <- this.ccm.xy$LibMeans$LibSize
    # Compute mean and sd rho from the largest MEAN.FRAC libraries
    min.libsize  <- this.ccm.xy$LibMeans$LibSize[ ceiling( nrow( this.ccm.xy$LibMeans ) * (1-MEAN.FRAC) ) ]
    max.libsize  <- max( this.ccm.xy$LibMeans$LibSize )
    rho.libsizes <- this.ccm.xy$LibMeans$LibSize[ which( this.ccm.xy$LibMeans$LibSize >= min.libsize &
                                                           this.ccm.xy$LibMeans$LibSize <= max.libsize ) ] |> unique()
    rho.rows     <- which( this.ccm.xy$CCM1_PredictStat$LibSize %in% rho.libsizes )
    # Compute lib:target rho mean & sd
    this.rho.xy <- mean( this.ccm.xy$CCM1_PredictStat$rho[ rho.rows ] )
    rho.xy[ which( DELAY.SEQ==i ) ] <- this.rho.xy
    this.sd.xy <- sd( this.ccm.xy$CCM1_PredictStat$rho[ rho.rows ] )
    rho.sd.xy[ which( DELAY.SEQ==i ) ] <- this.sd.xy
    # Compute target:lib rho mean & sd
    this.rho.yx <- mean( this.ccm.yx$CCM2_PredictStat$rho[ rho.rows ] )
    rho.yx[ which( DELAY.SEQ==i ) ] <- this.rho.yx
    this.sd.yx <- sd( this.ccm.yx$CCM2_PredictStat$rho[ rho.rows ] )
    rho.sd.yx[ which( DELAY.SEQ==i ) ] <- this.sd.yx
    
    if( i > 0 ){ next }  # If current delay is positive, go to next iteration to skip plotting rho vs. libsize
    # Plot for rho vs. library size
    # Plot points for all CCM tests
    plot( x = this.ccm.xy$CCM1_PredictStat$LibSize,  # lib:target
          y = this.ccm.xy$CCM1_PredictStat$rho,
          pch=16, col=rgb(1,0.2,0,0.1), cex=0.5, ylim=c(-0.5,1),
          main = paste(var.response,'&',var.driver,'\nLag =',i),
          xlab = 'Library size', ylab = 'rho' )
    points( x = this.ccm.yx$CCM2_PredictStat$LibSize,  # target:lib
            y = this.ccm.yx$CCM2_PredictStat$rho,
            pch=17, col=rgb(0,0.3,1,0.1), cex=0.5 )
    # Plot mean rho curves
    lines( x = lib.seq, y = this.ccm.xy$LibMeans$`lib:target`, lwd = 3, col = rgb(1,0.2,0,0.6) )  # lib:target
    lines( x = lib.seq, y = this.ccm.yx$LibMeans$`target:lib`, lwd = 3, col = rgb(0,0.3,1,0.6) )   # target:lib
    # Points
    # points( x = lib.seq, y = this.ccm.xy$LibMeans$`lib:target`, lwd = 3, col = rgb(1,0.2,0,0.8), pch=16 ) # lib:target
    # points( x = lib.seq, y = this.ccm.yx$LibMeans$`target:lib`, lwd = 3, col = rgb(0,0,0,0.5), pch=16 ) # target:lib
    # Lines and text
    abline( h=seq(-1,1,0.1), col=rgb(0,0,0,0.2) )  # horizontal lines
    abline( h=0, col=1, lwd=1 ) # zero line
    abline( h=abs(corr[which(DELAY.SEQ==i)]), lty=3, lwd=2,  # correlation absolute value 
            col=if(corr[which(DELAY.SEQ==i)]>0){rgb(0,0,0,0.7)}else{rgb(0.5,0,0,0.6)} )
    legend( 'bottom', bty='n', text.font=2, text.col = c( rgb(1,0.2,0,0.8), rgb(0,0.3,1,0.8) ),
            legend = c( paste0( var.response,' xmap ',var.driver,': rho = ',round(rho.xy[which(DELAY.SEQ==i)],3)),
                        paste0( var.driver,' xmap ',var.response,': rho = ',round(rho.yx[which(DELAY.SEQ==i)],3)) )
    )
    
  }  # // end i loop(DELAY.SEQ)
  mtext( paste(dates[1],'\u2013',dates[2]), side = 1, line = -1, outer = TRUE )
  if(write.png){ dev.off() }
  
  # Plot for mean (sd) rho vs. delay
  if(write.png){ 
    png( paste0(gpath,'ccm__',var.response,'_X_',var.driver,'__',dates[1],'--',dates[2],'__rho-vs-delay.png'),
         height=5, width=5, res=600, units='in' ) }  # // end if(write.png)
  par(mfrow=c(1,1))
  # Plot mean rho for lib:target
  plot( x=DELAY.SEQ, y=rho.xy, ylim=c(-0.5,1), lwd=0.5, cex=0.5, pch=16, col=rgb(1,0.2,0,0.8),
        main = paste(var.response,'&',var.driver,'\n', dates[1], "\u2013", dates[2] ),
        xlab='Lag', ylab='rho' )
  lines(  x=DELAY.SEQ, y=rho.xy, lwd=2, col=rgb(1,0.2,0,1) )
  # Plot mean rho for target:lib
  points( x=DELAY.SEQ, y=rho.yx, ylim=c(0,1), lwd=0.5, cex=0.5, pch=16, col=rgb(0,0.3,1,0.8) )
  lines(  x=DELAY.SEQ, y=rho.yx, lwd=2, col=rgb(0,0.3,1,1) )
  # Plot sd bars
  for( i in 1:length(DELAY.SEQ) ){
    # lib:target
    arrows( x0 = DELAY.SEQ[i],
            y0 = rho.xy[i]+rho.sd.xy[i], y1 = rho.xy[i]-rho.sd.xy[i],
            lwd = 2, code = 3, length = 0.05, angle = 90, col=rgb(1,0.2,0,0.8) )
    # target:lib
    arrows( x0 = DELAY.SEQ[i],
            y0 = rho.yx[i]+rho.sd.yx[i], y1 = rho.yx[i]-rho.sd.yx[i],
            lwd = 2, code = 3, length = 0.05, angle = 90, col=rgb(0,0.3,1,0.8) )
  } #  end i loop(DELAY.SEQ)
  # Lines and text
  abline( v=0, lty=2, lwd=2, col=rgb(0,0,0,0.4) )  # vertical line at delay=0
  abline( h=seq(-1,1,0.1), col=rgb(0,0,0,0.2) )  # horizontal lines
  abline( h=0, col=1, lwd=2 ) # zero line
  legend( 'bottomright', bty='n', text.font=2,
          text.col = c( rgb(1,0.2,0,1),rgb(0,0.3,1,0.8) ),
          legend = c( paste0( var.response,' xmap ',var.driver ),
                      paste0( var.driver,' xmap ',var.response ) )
  )
  if(write.png){ dev.off() }
  
  # Format output
  out.xy <- data.frame( delay = DELAY.SEQ,
                        rho.mean.xy = rho.xy, 
                        rho.sd.xy = rho.sd.xy )
  out.yx <- data.frame( delay = DELAY.SEQ,
                        rho.mean.yx = rho.yx, 
                        rho.sd.yx = rho.sd.yx )
  
  return( list( rho.xy = out.xy, rho.yx = out.yx ) )
  
}  # // end plot.CCM()



#  Define function to test significance of rho values returned by plot.CCM()
##  Most arguments default to parameters specified at the top of this script. Arguments:
##      CCM.OUTPUT  list object containing delay.CCM() output
##      CCM.RHO     list object containing plot.CCM() output
##      DELAY       delay whose corresponding rho value will be tested for significance (integer)
##      LIB.STEP    library step size (integer)
##      MEAN.FRAC   the fraction of libraries from which to compute mean rho values of surrogates (numeric)
##      NSURR       the number of surrogates to generate for each test (integer)
##      METHOD      surrogate generation method (see SurrogateData() documentation for options)
##      PERIOD      Seasonality period for METHOD='seasonal' surrogates (T_period argument)
##      SEED        random seed
pval.CCM <- function( CCM.OUTPUT, CCM.RHO, DELAY, LIB.STEP=lib.step, MEAN.FRAC=ccm.frac,
                      NSURR=nsurr, METHOD='seasonal', PERIOD=12, SEED=2930174  ){
  
  this.data <- CCM.OUTPUT$data
  this.rho  <- CCM.RHO$rho.xy$rho.mean.xy[ which(CCM.RHO$rho.xy$delay==DELAY) ]
  d         <- CCM.OUTPUT$d     # embedding delay (from library)
  tw        <- CCM.OUTPUT$tw   # Theiler window (from library)
  m         <- CCM.OUTPUT$m   # embedding dimension (from library)
  
  # Generate surrogates (each surrogate vector is a column in surr.data matrix)
  surr.mat <- SurrogateData( ts = this.data[,2],
                             method = METHOD, T_period = PERIOD,
                             num_surr = NSURR )
  
  # Define function to run CCM using a single surrogate library
  surr.CCM <- function( SURR, DF=this.data, FRAC=MEAN.FRAC ){
    
    # Build input dataframe 
    input.data <- DF  # library and target vectors from CCM.OUTPUT
    input.data$lib <- SURR  # replace library vector with surrogate
    # Adjust library and target vectors for delay
    if( DELAY<0 ){  # Negative delay: Response (library) predicts past values of driver (target)
      this.lib    <- input.data[ (-DELAY+1):nrow(input.data) , 2 ] # for delay=-1, lib takes 2nd row (-(-1)+1=2) thru final row..
      this.target <- input.data[ 1:(nrow(input.data)+DELAY)  , 3 ] # ..and target takes 1st row thru penultimate row (nrow+(-1)).
      surr.block   <- data.frame( t=1:length(this.lib), lib=this.lib, target=this.target )
    } else if ( DELAY>0 ){  # Positive delay: Response (library) predicts future values of driver (target)
      this.lib    <- input.data[ 1:(nrow(input.data)-DELAY) , 2 ] # for delay=1, lib takes first row through penultimate row..
      this.target <- input.data[ (DELAY+1):nrow(input.data) , 3 ] # ..and target takes 2nd row (1+1=2) through final row.
      surr.block   <- data.frame( t=1:length(this.lib), lib = this.lib, target = this.target )
    } else {  # No delay: Response predicts concurrent values of driver
      surr.block <- input.data
    }  # // end if(delay)
    
    # Run CCM to get surrogate rho values
    lib.max <- nrow(surr.block) - (m-1)*d  # max library size is constrained to this value when replacement=FALSE in CCM() call
    lib.seq <- seq( lib.max%%LIB.STEP, lib.max, LIB.STEP )
    if( lib.seq[1]<LIB.STEP ){ lib.seq <- lib.seq[-1] }
    # Run CCM using i'th surrogate as library
    this.CCM <-  CCM( dataFrame = surr.block, E = m, tau = -d, Tp = 0, exclusionRadius = tw,
                      columns = 'lib', target = 'target',
                      libSizes = lib.seq, sample = 1,
                      seed = SEED, includeData = FALSE, showPlot = FALSE )
    
    # Compute surrogate convergence value as the mean rho among the largest FRAC of libraries
    rho.rows <- ceiling( length( this.CCM$LibSize ) * (1-FRAC) ) : length(this.CCM$LibSize)
    this.mean.rho <- mean( this.CCM$`lib:target`[ rho.rows ] )
    
    return( this.mean.rho )
    
  }  # // end surr.CCM()
  
  # Run CCM using each surrogate vector as library
  surr.out <- apply( surr.mat, 2, surr.CCM )  # surr.out stores rho values from each surrogate test
  
  # Compute p-value
  k <- which( surr.out >= this.rho ) |> length()  # number of surrogate rho values exceeding the library's rho
  pval <- (k+1)/(NSURR+1)  # p-value: (k+1)/(n+1) where k is number of 'successes'
  
  return( list( pval = pval, ccm.rho = this.rho, surr.rho = surr.out, delay = DELAY ) )
  
}  # // end pval.CCM()


## Define function to generate plots
plot.CCM.results <- function( CCM.OUTPUT, CCM.RHO, DELAY.SEQ=delay.seq, PVAL, MAIN,
                              XLAB = FALSE, YLAB = FALSE, YMIN = 0.75, XSEQ = c(-7:0) ){
  
  # Get DELAY.SEQ indices
  delay.idx <- which(DELAY.SEQ %in% XSEQ)
  
  # Mean rho vs. lag
  plot( x=DELAY.SEQ[delay.idx],
        y=CCM.RHO$rho.xy$rho.mean.xy[delay.idx],
        xlab = "", ylab = "",ylim = c( YMIN,1 ),
        col=rgb(0,0,0,0), las = 1, cex.axis = 1.3 )
  mtext( MAIN, side = 3, line = 0, font = 1, cex = 1.3, adj = 0 )
  if(YLAB){ mtext( "rho", side = 2, cex = 0.9, line = 3.5 ) }
  if(XLAB){ mtext( "Delay (months)", side = 1, cex = 0.9, line = 2 ) }
  lines(  x=DELAY.SEQ[delay.idx],
          y=CCM.RHO$rho.xy$rho.mean.xy[delay.idx],
          lwd=2, col=rgb(1,0.2,0,0.7) )
  
  # st. dev. bars
  for( i in 1:length(DELAY.SEQ[delay.idx]) ){
    # lib:target
    if(CCM.RHO$rho.xy$rho.sd.xy[i] != 0 ){
      arrows( x0 = DELAY.SEQ[i],
              y0 = CCM.RHO$rho.xy$rho.mean.xy[i]+CCM.RHO$rho.xy$rho.sd.xy[i],
              y1 = CCM.RHO$rho.xy$rho.mean.xy[i]-CCM.RHO$rho.xy$rho.sd.xy[i],
              lwd = 0.5, code = 3, length = 0.03, angle = 90, col=rgb(0,0,0,1) )
    } else {
      points( x = DELAY.SEQ[i], y = CCM.RHO$rho.xy$rho.mean.xy[i],
              pch = "_" )
    }
  } #  end i loop(DELAY.SEQ)
  # Gridlines
  if( max(DELAY.SEQ[delay.idx])>0 ){
    abline( v=0, lty=2, lwd=1.5, col=rgb(0,0,0,0.4) )  # vertical line at delay=0
  }
  # abline( h=seq(-1,1,0.01), col=rgb(0,0,0,0.1) )  # horizontal lines
  # abline( h=seq(-1,1,0.1), col=rgb(0,0,0,0.1) )  # horizontal lines
  abline( h=axTicks(2), col=rgb(0,0,0,0.1) ) # horizontal lines
  abline( h=seq(min(axTicks(2)),max(axTicks(2)),diff(axTicks(2)[1:2])/2),
          col=rgb(0,0,0,0.1) ) # horizontal lines
  abline( h=0, col=1, lwd=2 ) # zero line
  
  # Significance (p-values from surrogate tests)
  sig <- which( PVAL <=0.05 )
  points( x=seq(0,min(DELAY.SEQ),-1)[sig], y=rep(1,length(PVAL))[sig],
          pch='*', lwd=1, cex=2.5, col=rgb(1,0.2,0,0.8) )
  
}  # // end plot.CCM.results()


# CCM tests ---------------------------------------------------------------


# OTB chl xmap Pyro
ccm.otb.chl.pyro <- delay.CCM( LIBRARY='chl_otb', TARGET='pyro_otb',
                             EMB = as.list(embed_tbl[which(embed_tbl$var=="chl_otb"),]) )
ccm.otb.chl.pyro.rho <- plot.CCM( ccm.otb.chl.pyro )
ccm.otb.chl.pyro.pval.d0  <- pval.CCM( CCM.OUTPUT=ccm.otb.chl.pyro, CCM.RHO=ccm.otb.chl.pyro.rho, DELAY = 0 )
ccm.otb.chl.pyro.pval.d1  <- pval.CCM( CCM.OUTPUT=ccm.otb.chl.pyro, CCM.RHO=ccm.otb.chl.pyro.rho, DELAY = -1 ) 
ccm.otb.chl.pyro.pval.d2  <- pval.CCM( CCM.OUTPUT=ccm.otb.chl.pyro, CCM.RHO=ccm.otb.chl.pyro.rho, DELAY = -2 ) 
ccm.otb.chl.pyro.pval.d3  <- pval.CCM( CCM.OUTPUT=ccm.otb.chl.pyro, CCM.RHO=ccm.otb.chl.pyro.rho, DELAY = -3 ) 
ccm.otb.chl.pyro.pval.d4  <- pval.CCM( CCM.OUTPUT=ccm.otb.chl.pyro, CCM.RHO=ccm.otb.chl.pyro.rho, DELAY = -4 ) 
ccm.otb.chl.pyro.pval.d5  <- pval.CCM( CCM.OUTPUT=ccm.otb.chl.pyro, CCM.RHO=ccm.otb.chl.pyro.rho, DELAY = -5 ) 
ccm.otb.chl.pyro.pval.d6  <- pval.CCM( CCM.OUTPUT=ccm.otb.chl.pyro, CCM.RHO=ccm.otb.chl.pyro.rho, DELAY = -6 ) 
ccm.otb.chl.pyro.pval.d7  <- pval.CCM( CCM.OUTPUT=ccm.otb.chl.pyro, CCM.RHO=ccm.otb.chl.pyro.rho, DELAY = -7 ) 

# NW chl xmap Pyro
ccm.nw.chl.pyro <- delay.CCM( LIBRARY='chl_otb_nw', TARGET='pyro_otb_nw',
                            EMB = as.list(embed_tbl[which(embed_tbl$var=="chl_otb_nw"),]) )
ccm.nw.chl.pyro.rho <- plot.CCM( ccm.nw.chl.pyro )
ccm.nw.chl.pyro.pval.d0  <- pval.CCM( CCM.OUTPUT=ccm.nw.chl.pyro, CCM.RHO=ccm.nw.chl.pyro.rho, DELAY = 0 )
ccm.nw.chl.pyro.pval.d1  <- pval.CCM( CCM.OUTPUT=ccm.nw.chl.pyro, CCM.RHO=ccm.nw.chl.pyro.rho, DELAY = -1 ) 
ccm.nw.chl.pyro.pval.d2  <- pval.CCM( CCM.OUTPUT=ccm.nw.chl.pyro, CCM.RHO=ccm.nw.chl.pyro.rho, DELAY = -2 ) 
ccm.nw.chl.pyro.pval.d3  <- pval.CCM( CCM.OUTPUT=ccm.nw.chl.pyro, CCM.RHO=ccm.nw.chl.pyro.rho, DELAY = -3 ) 
ccm.nw.chl.pyro.pval.d4  <- pval.CCM( CCM.OUTPUT=ccm.nw.chl.pyro, CCM.RHO=ccm.nw.chl.pyro.rho, DELAY = -4 ) 
ccm.nw.chl.pyro.pval.d5  <- pval.CCM( CCM.OUTPUT=ccm.nw.chl.pyro, CCM.RHO=ccm.nw.chl.pyro.rho, DELAY = -5 ) 
ccm.nw.chl.pyro.pval.d6  <- pval.CCM( CCM.OUTPUT=ccm.nw.chl.pyro, CCM.RHO=ccm.nw.chl.pyro.rho, DELAY = -6 ) 
ccm.nw.chl.pyro.pval.d7  <- pval.CCM( CCM.OUTPUT=ccm.nw.chl.pyro, CCM.RHO=ccm.nw.chl.pyro.rho, DELAY = -7 ) 

# NE chl xmap Pyro
ccm.ne.chl.pyro <- delay.CCM( LIBRARY='chl_otb_ne', TARGET='pyro_otb_ne',
                            EMB = as.list(embed_tbl[which(embed_tbl$var=="chl_otb_ne"),]) )
ccm.ne.chl.pyro.rho <- plot.CCM( ccm.ne.chl.pyro )
ccm.ne.chl.pyro.pval.d0  <- pval.CCM( CCM.OUTPUT=ccm.ne.chl.pyro, CCM.RHO=ccm.ne.chl.pyro.rho, DELAY = 0 )
ccm.ne.chl.pyro.pval.d1  <- pval.CCM( CCM.OUTPUT=ccm.ne.chl.pyro, CCM.RHO=ccm.ne.chl.pyro.rho, DELAY = -1 ) 
ccm.ne.chl.pyro.pval.d2  <- pval.CCM( CCM.OUTPUT=ccm.ne.chl.pyro, CCM.RHO=ccm.ne.chl.pyro.rho, DELAY = -2 ) 
ccm.ne.chl.pyro.pval.d3  <- pval.CCM( CCM.OUTPUT=ccm.ne.chl.pyro, CCM.RHO=ccm.ne.chl.pyro.rho, DELAY = -3 ) 
ccm.ne.chl.pyro.pval.d4  <- pval.CCM( CCM.OUTPUT=ccm.ne.chl.pyro, CCM.RHO=ccm.ne.chl.pyro.rho, DELAY = -4 ) 
ccm.ne.chl.pyro.pval.d5  <- pval.CCM( CCM.OUTPUT=ccm.ne.chl.pyro, CCM.RHO=ccm.ne.chl.pyro.rho, DELAY = -5 ) 
ccm.ne.chl.pyro.pval.d6  <- pval.CCM( CCM.OUTPUT=ccm.ne.chl.pyro, CCM.RHO=ccm.ne.chl.pyro.rho, DELAY = -6 ) 
ccm.ne.chl.pyro.pval.d7  <- pval.CCM( CCM.OUTPUT=ccm.ne.chl.pyro, CCM.RHO=ccm.ne.chl.pyro.rho, DELAY = -7 )

# CW chl xmap Pyro
ccm.cw.chl.pyro <- delay.CCM( LIBRARY='chl_otb_cw', TARGET='pyro_otb_cw',
                            EMB = as.list(embed_tbl[which(embed_tbl$var=="chl_otb_cw"),]) )
ccm.cw.chl.pyro.rho <- plot.CCM( ccm.cw.chl.pyro )
ccm.cw.chl.pyro.pval.d0  <- pval.CCM( CCM.OUTPUT=ccm.cw.chl.pyro, CCM.RHO=ccm.cw.chl.pyro.rho, DELAY = 0 )
ccm.cw.chl.pyro.pval.d1  <- pval.CCM( CCM.OUTPUT=ccm.cw.chl.pyro, CCM.RHO=ccm.cw.chl.pyro.rho, DELAY = -1 ) 
ccm.cw.chl.pyro.pval.d2  <- pval.CCM( CCM.OUTPUT=ccm.cw.chl.pyro, CCM.RHO=ccm.cw.chl.pyro.rho, DELAY = -2 ) 
ccm.cw.chl.pyro.pval.d3  <- pval.CCM( CCM.OUTPUT=ccm.cw.chl.pyro, CCM.RHO=ccm.cw.chl.pyro.rho, DELAY = -3 ) 
ccm.cw.chl.pyro.pval.d4  <- pval.CCM( CCM.OUTPUT=ccm.cw.chl.pyro, CCM.RHO=ccm.cw.chl.pyro.rho, DELAY = -4 ) 
ccm.cw.chl.pyro.pval.d5  <- pval.CCM( CCM.OUTPUT=ccm.cw.chl.pyro, CCM.RHO=ccm.cw.chl.pyro.rho, DELAY = -5 ) 
ccm.cw.chl.pyro.pval.d6  <- pval.CCM( CCM.OUTPUT=ccm.cw.chl.pyro, CCM.RHO=ccm.cw.chl.pyro.rho, DELAY = -6 ) 
ccm.cw.chl.pyro.pval.d7  <- pval.CCM( CCM.OUTPUT=ccm.cw.chl.pyro, CCM.RHO=ccm.cw.chl.pyro.rho, DELAY = -7 ) 

# CE chl xmap Pyro
ccm.ce.chl.pyro <- delay.CCM( LIBRARY='chl_otb_ce', TARGET='pyro_otb_ce',
                            EMB = as.list(embed_tbl[which(embed_tbl$var=="chl_otb_ce"),]) )
ccm.ce.chl.pyro.rho <- plot.CCM( ccm.ce.chl.pyro )
ccm.ce.chl.pyro.pval.d0  <- pval.CCM( CCM.OUTPUT=ccm.ce.chl.pyro, CCM.RHO=ccm.ce.chl.pyro.rho, DELAY = 0 )
ccm.ce.chl.pyro.pval.d1  <- pval.CCM( CCM.OUTPUT=ccm.ce.chl.pyro, CCM.RHO=ccm.ce.chl.pyro.rho, DELAY = -1 ) 
ccm.ce.chl.pyro.pval.d2  <- pval.CCM( CCM.OUTPUT=ccm.ce.chl.pyro, CCM.RHO=ccm.ce.chl.pyro.rho, DELAY = -2 ) 
ccm.ce.chl.pyro.pval.d3  <- pval.CCM( CCM.OUTPUT=ccm.ce.chl.pyro, CCM.RHO=ccm.ce.chl.pyro.rho, DELAY = -3 ) 
ccm.ce.chl.pyro.pval.d4  <- pval.CCM( CCM.OUTPUT=ccm.ce.chl.pyro, CCM.RHO=ccm.ce.chl.pyro.rho, DELAY = -4 ) 
ccm.ce.chl.pyro.pval.d5  <- pval.CCM( CCM.OUTPUT=ccm.ce.chl.pyro, CCM.RHO=ccm.ce.chl.pyro.rho, DELAY = -5 ) 
ccm.ce.chl.pyro.pval.d6  <- pval.CCM( CCM.OUTPUT=ccm.ce.chl.pyro, CCM.RHO=ccm.ce.chl.pyro.rho, DELAY = -6 ) 
ccm.ce.chl.pyro.pval.d7  <- pval.CCM( CCM.OUTPUT=ccm.ce.chl.pyro, CCM.RHO=ccm.ce.chl.pyro.rho, DELAY = -7 ) 

# SW chl xmap Pyro
ccm.sw.chl.pyro <- delay.CCM( LIBRARY='chl_otb_sw', TARGET='pyro_otb_sw',
                            EMB = as.list(embed_tbl[which(embed_tbl$var=="chl_otb_sw"),]) )
ccm.sw.chl.pyro.rho <- plot.CCM( ccm.sw.chl.pyro )
ccm.sw.chl.pyro.pval.d0  <- pval.CCM( CCM.OUTPUT=ccm.sw.chl.pyro, CCM.RHO=ccm.sw.chl.pyro.rho, DELAY = 0 )
ccm.sw.chl.pyro.pval.d1  <- pval.CCM( CCM.OUTPUT=ccm.sw.chl.pyro, CCM.RHO=ccm.sw.chl.pyro.rho, DELAY = -1 ) 
ccm.sw.chl.pyro.pval.d2  <- pval.CCM( CCM.OUTPUT=ccm.sw.chl.pyro, CCM.RHO=ccm.sw.chl.pyro.rho, DELAY = -2 ) 
ccm.sw.chl.pyro.pval.d3  <- pval.CCM( CCM.OUTPUT=ccm.sw.chl.pyro, CCM.RHO=ccm.sw.chl.pyro.rho, DELAY = -3 ) 
ccm.sw.chl.pyro.pval.d4  <- pval.CCM( CCM.OUTPUT=ccm.sw.chl.pyro, CCM.RHO=ccm.sw.chl.pyro.rho, DELAY = -4 ) 
ccm.sw.chl.pyro.pval.d5  <- pval.CCM( CCM.OUTPUT=ccm.sw.chl.pyro, CCM.RHO=ccm.sw.chl.pyro.rho, DELAY = -5 ) 
ccm.sw.chl.pyro.pval.d6  <- pval.CCM( CCM.OUTPUT=ccm.sw.chl.pyro, CCM.RHO=ccm.sw.chl.pyro.rho, DELAY = -6 ) 
ccm.sw.chl.pyro.pval.d7  <- pval.CCM( CCM.OUTPUT=ccm.sw.chl.pyro, CCM.RHO=ccm.sw.chl.pyro.rho, DELAY = -7 ) 

# SE chl xmap Pyro
ccm.se.chl.pyro <- delay.CCM( LIBRARY='chl_otb_se', TARGET='pyro_otb_se',
                            EMB = as.list(embed_tbl[which(embed_tbl$var=="chl_otb_se"),]) )
ccm.se.chl.pyro.rho <- plot.CCM( ccm.se.chl.pyro )
ccm.se.chl.pyro.pval.d0  <- pval.CCM( CCM.OUTPUT=ccm.se.chl.pyro, CCM.RHO=ccm.se.chl.pyro.rho, DELAY = 0 )
ccm.se.chl.pyro.pval.d1  <- pval.CCM( CCM.OUTPUT=ccm.se.chl.pyro, CCM.RHO=ccm.se.chl.pyro.rho, DELAY = -1 ) 
ccm.se.chl.pyro.pval.d2  <- pval.CCM( CCM.OUTPUT=ccm.se.chl.pyro, CCM.RHO=ccm.se.chl.pyro.rho, DELAY = -2 ) 
ccm.se.chl.pyro.pval.d3  <- pval.CCM( CCM.OUTPUT=ccm.se.chl.pyro, CCM.RHO=ccm.se.chl.pyro.rho, DELAY = -3 ) 
ccm.se.chl.pyro.pval.d4  <- pval.CCM( CCM.OUTPUT=ccm.se.chl.pyro, CCM.RHO=ccm.se.chl.pyro.rho, DELAY = -4 ) 
ccm.se.chl.pyro.pval.d5  <- pval.CCM( CCM.OUTPUT=ccm.se.chl.pyro, CCM.RHO=ccm.se.chl.pyro.rho, DELAY = -5 ) 
ccm.se.chl.pyro.pval.d6  <- pval.CCM( CCM.OUTPUT=ccm.se.chl.pyro, CCM.RHO=ccm.se.chl.pyro.rho, DELAY = -6 ) 
ccm.se.chl.pyro.pval.d7  <- pval.CCM( CCM.OUTPUT=ccm.se.chl.pyro, CCM.RHO=ccm.se.chl.pyro.rho, DELAY = -7 ) 



# CCM results -------------------------------------------------------------

# Report rho values and p-values
rhovars <- paste0( ls()[ grep('pval.d',ls()) ], '$ccm.rho' ) |> data.frame() # find all rho variable names in env
rhovals <- apply( rhovars, 1, function(x) eval(parse(text=x)) )  # get all rho values
pvars <- paste0( ls()[ grep('pval.d',ls()) ], '$pval' ) |> data.frame()  # find all pval variable names in env
pvals <- apply( pvars, 1, function(x) eval(parse(text=x)) )  # get all p-values
# Define function to assign stars to various significance levels
sig.stars <- function(X,a=c(0.10,0.05,0.01,0.001)){
  if(X<a[4]){  # <0.001
    stars <- '***'
  } else if(X<a[3]){  # <0.01
    stars <- '**'
  } else if(X<a[2]){  # <0.05
    stars <- '*'
  } else if(X<a[1]){  #  <0.10
    stars <- '.'
  } else {
    stars <- ''
  }
  return(stars)
}  # // end sig.stars()
# Construct dataframe to store p-values
df.pval <- data.frame( var1 = sub('[\\$,]pval','',as.vector(unlist(pvars))),
                       rho = round(rhovals,3),
                       pval = round(pvals,3),
                       sig  = apply( data.frame(pvals), 1, sig.stars )
)
df.pval[ which(df.pval$pval<0.1), ]  # Print 'significant' results to console

# Export outputs
save.image( file = "../data/ccm_chl_pyro.RData" )
write.csv( df.pval, file = "../data/ccm_chl_pyro_pvals.csv", row.names = FALSE )


# CCM figures --------------------------------------------------------------
png( "../figs/ccm_chl_pyro.png", width = 6, height = 10, units = 'in', res = 600 )
par(mfrow=c(4,2), mar=c(4,4.5,2,1))

plot.CCM.results( CCM.OUTPUT=ccm.otb.chl.pyro,
                  CCM.RHO=ccm.otb.chl.pyro.rho,
                  MAIN = "OTB segment", YLAB = TRUE,
                  PVAL=c(ccm.otb.chl.pyro.pval.d0$pval, 
                         ccm.otb.chl.pyro.pval.d1$pval, 
                         ccm.otb.chl.pyro.pval.d2$pval, 
                         ccm.otb.chl.pyro.pval.d3$pval, 
                         ccm.otb.chl.pyro.pval.d4$pval, 
                         ccm.otb.chl.pyro.pval.d5$pval, 
                         ccm.otb.chl.pyro.pval.d6$pval, 
                         ccm.otb.chl.pyro.pval.d7$pval ) )
plot(c(0,1),c(0,1),bty='n',xaxt='n',yaxt='n',xlab='',ylab='',col=rgb(0,0,0,0))
text( x = 0, y = 0.95, pos = 4, cex = 1.3,
      labels = "Chlorophyll-a xmap Pyro" )
text( x = 0, y = 0.85, pos = 4, cex = 1.3,
      labels = paste0("(2012\U2013","2021)") )
segments( x0 = 0, x1 = 0.05, y0 = 0.6, lwd=2, col=rgb(1,0.2,0,0.7) )
text( x = 0.09, y = 0.6, pos = 4, labels = "xmap correlation\n(mean rho \U00B1SD)", cex = 1.3 )
arrows( x0 = 0.025, y0 = 0.57, y1 = 0.63, 
        lwd = 0.5, code = 3, length = 0.03, angle = 90 )
points( x = 0.025, y = 0.3, pch = "*", cex = 2.5, col=rgb(1,0.2,0,0.8) )
text( x = 0.09, y = 0.3, pos = 4, labels = "P<0.05\n(surrogate data test)", cex = 1.3 )

plot.CCM.results( CCM.OUTPUT=ccm.nw.chl.pyro,
                  CCM.RHO=ccm.nw.chl.pyro.rho,
                  MAIN = "NW sub-segment", YLAB = TRUE,
                  PVAL=c(ccm.nw.chl.pyro.pval.d0$pval, 
                         ccm.nw.chl.pyro.pval.d1$pval, 
                         ccm.nw.chl.pyro.pval.d2$pval, 
                         ccm.nw.chl.pyro.pval.d3$pval, 
                         ccm.nw.chl.pyro.pval.d4$pval, 
                         ccm.nw.chl.pyro.pval.d5$pval, 
                         ccm.nw.chl.pyro.pval.d6$pval, 
                         ccm.nw.chl.pyro.pval.d7$pval ) )
plot.CCM.results( CCM.OUTPUT=ccm.ne.chl.pyro,
                  CCM.RHO=ccm.ne.chl.pyro.rho,
                  MAIN = "NE sub-segment",
                  PVAL=c(ccm.ne.chl.pyro.pval.d0$pval, 
                         ccm.ne.chl.pyro.pval.d1$pval, 
                         ccm.ne.chl.pyro.pval.d2$pval, 
                         ccm.ne.chl.pyro.pval.d3$pval, 
                         ccm.ne.chl.pyro.pval.d4$pval, 
                         ccm.ne.chl.pyro.pval.d5$pval, 
                         ccm.ne.chl.pyro.pval.d6$pval, 
                         ccm.ne.chl.pyro.pval.d7$pval ) )
plot.CCM.results( CCM.OUTPUT=ccm.cw.chl.pyro,
                  CCM.RHO=ccm.cw.chl.pyro.rho,
                  MAIN = "CW sub-segment", YLAB = TRUE,
                  PVAL=c(ccm.cw.chl.pyro.pval.d0$pval, 
                         ccm.cw.chl.pyro.pval.d1$pval, 
                         ccm.cw.chl.pyro.pval.d2$pval, 
                         ccm.cw.chl.pyro.pval.d3$pval, 
                         ccm.cw.chl.pyro.pval.d4$pval, 
                         ccm.cw.chl.pyro.pval.d5$pval, 
                         ccm.cw.chl.pyro.pval.d6$pval, 
                         ccm.cw.chl.pyro.pval.d7$pval ) )
plot.CCM.results( CCM.OUTPUT=ccm.ce.chl.pyro,
                  CCM.RHO=ccm.ce.chl.pyro.rho,
                  MAIN = "CE sub-segment",
                  PVAL=c(ccm.ce.chl.pyro.pval.d0$pval, 
                         ccm.ce.chl.pyro.pval.d1$pval, 
                         ccm.ce.chl.pyro.pval.d2$pval, 
                         ccm.ce.chl.pyro.pval.d3$pval, 
                         ccm.ce.chl.pyro.pval.d4$pval, 
                         ccm.ce.chl.pyro.pval.d5$pval, 
                         ccm.ce.chl.pyro.pval.d6$pval, 
                         ccm.ce.chl.pyro.pval.d7$pval ) )
plot.CCM.results( CCM.OUTPUT=ccm.sw.chl.pyro,
                  CCM.RHO=ccm.sw.chl.pyro.rho,
                  MAIN = "SW sub-segment", YLAB = TRUE, XLAB = TRUE,
                  PVAL=c(ccm.sw.chl.pyro.pval.d0$pval, 
                         ccm.sw.chl.pyro.pval.d1$pval, 
                         ccm.sw.chl.pyro.pval.d2$pval, 
                         ccm.sw.chl.pyro.pval.d3$pval, 
                         ccm.sw.chl.pyro.pval.d4$pval, 
                         ccm.sw.chl.pyro.pval.d5$pval, 
                         ccm.sw.chl.pyro.pval.d6$pval, 
                         ccm.sw.chl.pyro.pval.d7$pval ) )
plot.CCM.results( CCM.OUTPUT=ccm.se.chl.pyro,
                  CCM.RHO=ccm.se.chl.pyro.rho,
                  MAIN = "SE sub-segment", XLAB = TRUE,
                  PVAL=c(ccm.se.chl.pyro.pval.d0$pval, 
                         ccm.se.chl.pyro.pval.d1$pval, 
                         ccm.se.chl.pyro.pval.d2$pval, 
                         ccm.se.chl.pyro.pval.d3$pval, 
                         ccm.se.chl.pyro.pval.d4$pval, 
                         ccm.se.chl.pyro.pval.d5$pval, 
                         ccm.se.chl.pyro.pval.d6$pval, 
                         ccm.se.chl.pyro.pval.d7$pval ) )

dev.off()