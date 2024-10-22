

# 
# ############################# U-Pb modeling inputs  ############################# 
library(IsoplotR)
library(dplyr)

## read in equations from this file
source( "UPb_Constants_Functions_Libraries.R" )   # Read in constants and functions from the other file



## read in the data file and save the sample name
sample.name <- "Sample_2_R" ## set the sample name
data.raw <- read.csv( paste( sample.name, ".csv", sep = "")) # read in sample file 

#read data into isoplotr
# now do the isoplot calculations
data.concordia <- IsoplotR::read.data(  data.raw[ , 2:6], ierr = 2,
                                   method = 'U-Pb', format = 1 )
concordia <- IsoplotR::concordia( data.concordia, type = 1  )

## now run our modeling reduction for the dataset
# Define the sample name, and node spacing
node.spacing	<-2

#############################  SWITCHES #############################
## this should be "Y" to normalize the uncertainties to the median value
#    otherwise it doesn't do anything
normalize.uncertainty	 <- "N"  

## this should be "detrital" to weight against concordant analyses
#    otherwise it should be 'single' to not weight against concordant analyses
data.type	 <- "single"  

## If cut.data.by.ratios is Y it trim the input data by the cuts below
cut.data.by.ratios	<- "N"
## These are the start and ends, in ratio space, this cuts data out of the data file
startcut.r75        <- 0
endcut.r75          <- 20
startcut.r68        <- 0
endcut.r68          <- 0.8

# This zooms the plots into a certain age window
#  		Use this to either simply zoom in on a particular age, or 
#  		to zoom in and use a very tight node spacing to save computational time	
#       it doesn't perform the analysis outside of the age window defined below
zoom.analysis		<- "Y"
## These are the start and ends, only performs the reduction on certain nodes defined
#	by the ages below here
startcut.age.lower        <- 0 			## Age in Ma
endcut.age.lower          <- 50		  ## Age in Ma
startcut.age.upper        <- 800 		## Age in Ma
endcut.age.upper          <- 1900		  ## Age in Ma


#Start resampling step
subsample.iteration <- 2 #number of subsamples, we used 1000 which will take awhile
## use a small number like 2 to test, then increase

#Make blank list to save results 
resample.datapoints <- list() #store bootstrapped U-Pb data
resample.results.data <- list() # store results from discordance dating

#resample loop to generate bootstrapped samples
for(j in 1:subsample.iteration){
  data.test <- dplyr::slice_sample(.data=data.raw, n=nrow(data.raw), replace=TRUE)
  resample.datapoints[[j]] <- data.test
  # setwd(source.dir)
  source( "UPb_Reduction_Resample.R" )  ## do the reduction for each iteration
  resample.results.data[[j]] <- lowerdisc.sum.total
  print(j)
}

# now use the resampled data and plot it
#peak finding function
list.maxposition <- lapply(resample.results.data,function(df)which.max(df$normalized.sum.likelihood))
list.maxlikelihood <- lapply(resample.results.data,function(df)max(df$normalized.sum.likelihood))
data.max.age <- data.frame( max.position=unlist(list.maxposition),
                            max.likelihood=unlist(list.maxlikelihood))
data.max.age$age <- data.max.age$max.position * node.spacing + startcut.age.lower


## expand the list
resampled.results.expanded <-  bind_rows( resample.results.data, .id = "column_label")
resampled.results.expanded$run.number <- resampled.results.expanded$column_label

## reimport raw data and do the reduction
data.test <- data.raw
source( "UPb_Reduction_Resample.R" )  ## do the reduction

## plot the lower interecept curves for models and the real curve
ggplot(dplyr::bind_rows(resample.results.data, .id="data_frame"),
       aes(x=`Lower Intercept`/1e6, y=normalized.sum.likelihood, color=data_frame)) +
  geom_line( show.legend = FALSE, alpha = 0.2 ) +
  xlab( 'Intercept Age (Ma)' ) +
  ylab( 'Normalized Likelihood' ) +
  # geom_vline( xintercept = 43.5 ) +
  fte_theme_white() + 
  scale_colour_grey( start = 0.3, end = 0.31 ) +
  geom_line( data = lowerdisc.sum.total, color = "black" )



#histogram of max values
ggplot(data.max.age,aes(x=age)) + 
    geom_histogram( bins = 50 ) + 
    xlim( 0, 50 ) +
    xlab( 'Maximum Age (Ma)' ) +
    ylab( 'Number' ) +
    fte_theme_white() 

#scatter of max values
  ggplot(data.max.age,aes(x=age, y=max.likelihood)) + 
    geom_point() + 
    xlim( 0, 50 ) +
    xlab( 'Maximum Age (Ma)' ) +
    ylab( 'Likelihood' ) +
    fte_theme_white() 

 



          