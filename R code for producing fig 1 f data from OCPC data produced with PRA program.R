# diQ4p567b for a rerun diQ5p559p1b.R (16_08_27) derived from diQ4p558b.R (16_08_27) derived from diQ4p553p5a (16_08_17)
# diQ3p442 derived from diQ3p431b derived from diQ3p407b (16_02_12)
# R script used for analyzing Dicty growth
# Specifically, given two .pls files from "PRA pulse analyzer", one being a control, and the other being a test file, will produce 
# three control corrected plots as well as the integral of the pulse height histogram. Said integral will then be plotted against time
# to give a plot proportional to the actual cell growth plot.
# Made by: Carl Franck Research Group @ Cornell University 

# Clear Workspace
rm(list=ls(all=TRUE))

# Load necessary packages 
library(ggplot2)
library(scales)
library(grid)
library(Hmisc)


# Set full directory for where control file is and where test file(s) are, again make sure to change all '\'  to '/'
ControlDirectory <- "C:/Users/Carl/Desktop/diQ4p567a/control/"
TestDirectory <- "C:/Users/Carl/Desktop/diQ4p567a/test/"
OutPutDirectory <- "C:/Users/Carl/Desktop/diQ4p567a/output/"

# Loading Function
calculateandplotPulseData <- function(direct, pulsesC.data,HeightMax,HeightMin,WidthMax,WidthMin,alpha,beta) {
  
  # Sets working directory and loads all files 
  setwd(direct)
  files.list <- list.files()
  files.total <- length(files.list)
  # Variable to keep track of the area under the curve
  files.integrals <- rep(0, files.total)
  files.datas <- rep(0, files.total)
  # Iterates through all files
  HitMatrix <- matrix(rep(0,1500),nrow=500,ncol=1)
  
  #Make Control into Histogram
  BreaksW <- seq(0,20,by=.1)
  BreaksH <- seq(HeightMin*100/2147483647,HeightMax*100/2147483647,by=.1)
  HistInfoWidthC <- hist(pulsesC.data$pulseWidthMs,breaks=BreaksW)
  HistInfoHeightC <- hist(pulsesC.data$pulseHeightNorm,breaks=BreaksH)
  height.data <- data.frame(
    ControlHeight=HistInfoHeightC$counts
  )
  width.data <- data.frame(
    ControlWidth=HistInfoWidthC$counts
  )
  
  
  for (i in seq(1, files.total)) {
    
    # Variable to keep track of the count that matches the condition
    modifiedCount <- 1 
    
    # Show progress bar
    print ((i/files.total)*100)
    
    # Open file for reading in binary mode
    raw.file <- file(files.list[i], "rb")
    
    # Read file header
    raw.header <- readBin(raw.file, integer(), 1)
    raw.totalpulses <- readBin(raw.file, integer(), 1)
    raw.blank <- readBin(raw.file, integer(), 1)
    raw.blank2 <- readBin(raw.file, integer(), 1)
    raw.acqTime <- readBin(raw.file, integer(), 1)
    raw.headerRemaining.bin <- readBin(raw.file, raw(), 4100)
    
    # Create empty vectors for storing pulse information
    # This will work much faster than expanding them on the fly.
    c.beginTimes <- rep(0, raw.totalpulses)
    c.heights <- rep(0, raw.totalpulses)
    c.coincs <- rep(0, raw.totalpulses)
    c.widths <- rep(0, raw.totalpulses)
    
    
    # Iterates through all pulses in one file
    for (j in seq(1, raw.totalpulses)) {
      pulse.beginTime <- readBin(raw.file, integer(), n=1, size=8)
      pulse.width <- readBin(raw.file, integer(), n=1, size=2)
      pulse.coinc <- readBin(raw.file, integer(), n=1, size=2)
      pulse.height <- readBin(raw.file, integer(), n=1, size=4)
      
      # Checks if pulse meets criteria of width and height, if so records it
      if (pulse.width < WidthMax && pulse.width > WidthMin && 
          pulse.height > HeightMin && pulse.height<HeightMax){
        c.beginTimes[modifiedCount] <- pulse.beginTime
        c.heights[modifiedCount] <- pulse.height
        c.coincs[modifiedCount] <- pulse.coinc
        c.widths[modifiedCount] <- pulse.width
        modifiedCount <- modifiedCount + 1
      }
      
    }
    
    # Trim the array
    c.beginTimes <- c.beginTimes[1:modifiedCount]
    c.heights <- c.heights[1:modifiedCount]
    c.coincs <- c.coincs[1:modifiedCount]
    c.widths <- c.widths[1:modifiedCount]
    
    
    
    # Organize this data into a data frame
    pulses.data <- data.frame(
      pulseId = seq(1, modifiedCount),
      pulseBeginTime = c.beginTimes,
      pulseHeight = c.heights,
      pulseCoincidence = c.coincs,
      pulseWidth = c.widths
    )

    
    # Normalize the raw numbers
    pulses.data$pulseHeightNorm <- pulses.data$pulseHeight*100/2147483647
    pulses.data$pulseWidthMs <- pulses.data$pulseWidth/48
    # THis is divided by 48 since the sampling rate is 48000 times per second
    
    # Correcting signal for control
    HistInfoWidthT <- hist(pulses.data$pulseWidthMs,breaks=BreaksW)
    ControlFixedCountsW <- HistInfoWidthT$counts - HistInfoWidthC$counts
    ErrorW <- sqrt(HistInfoWidthT$counts + HistInfoWidthC$counts)/2
    
    HistInfoHeightT <- hist(pulses.data$pulseHeightNorm,breaks=BreaksH)
    ControlFixedCountsH <- HistInfoHeightT$counts - HistInfoHeightC$counts
    ErrorH <- sqrt(HistInfoHeightT$counts + HistInfoHeightC$counts)/2
    
    height.data[i+1]=HistInfoHeightT$counts
    width.data[i+1]=HistInfoWidthT$counts
    
    setwd(OutPutDirectory)
    
    files.integrals[i] <- sum(ControlFixedCountsH)
    #files.integrals[i] <- length(pulses.data$pulseHeightNorm)
    
    if (alpha == 1){
    #   Make plots
    # Width Histogram
    dir <- paste(OutPutDirectory,i,sep="")
    dir <- paste(dir,"W.jpeg",sep="")
    jpeg(file=dir)
    plot(HistInfoWidthT$mids,ControlFixedCountsW, xlim=c(0,30), ylim=c(-250,250), type="o")
    
    dev.off()
    
    # Height Histogram
    dir <- paste(OutPutDirectory,i,sep="")
    dir <- paste(dir,"H.jpeg",sep="")
    jpeg(file=dir)
    # changing height hist plot limits xlim=c(0,100) changed to (0,25) ylim=(-10,20) changed to (-5,40)
    # for late run, changing vertical range from ylim=c(-5.40) to c(-5,300) and horiz range from (0,25) to (0,75)
    plot(HistInfoHeightT$mids,ControlFixedCountsH, xlim=c(0,75), ylim=c(-5,300), type="o")
    dev.off()
    
    # Heat map of Width vs Height and counts
    dir <- paste(OutPutDirectory,i,sep="")
    dir <- paste(dir,"M.jpeg",sep="")
    jpeg(file=dir)
    #   
    p <- ggplot(pulses.data, aes(x=pulseHeightNorm, y=pulseWidthMs))
    p + geom_point()
    #   
    print(
  #adjusting the heatmap intensity scale
  # changed pulse width range from (0,30) to (0,.5), changed ratio from .5 to 20
      p + theme_bw() + stat_bin2d(drop=F, binwidth=c(0.04, 0.06)) + 
        scale_x_continuous("pulse height (a.u.)", expand=c(0,0), limits=c(0,50)) +
        scale_y_continuous("pulse width (ms)", expand=c(0,0), limits=c(0,.5)) + 
        coord_fixed(ratio=20) +
        scale_fill_gradient(
 # changed the na. value from #000000
          trans="log10", na.value="#ffffff", 
          low="#ffffff", high="#000000",
  #       low="#000000", high="#71ff71",
          limits=c(0.2,1)
        ) +
        theme(
          panel.grid = element_blank()   
        )
    )
    }
    Temp = sum(matrix(HistInfoHeightT$counts,nrow=500,ncol=1))
    TempM = matrix(HistInfoHeightT$counts,nrow=500,ncol=1)/Temp
    HitMatrix <- cbind(HitMatrix,TempM)
    
    dev.off()
    setwd(direct)
  }
  if (beta == 3){
  setwd(OutPutDirectory)
  #Results$integrals <- Results$integrals - alpha
  write.csv(width.data,"WidthHisto.csv")
  write.csv(height.data,"HeightHisto.csv")
  # Returns integral vector
  ListReturn <- list("integrals" = files.integrals, "HitMatrix" = HitMatrix, "Counts" = ControlFixedCountsH)
  return(ListReturn)
  }
}

calculateControl <- function(direct,HeightMax,HeightMin,WidthMax,WidthMin)
{ 
  # Sets working directory and loads all files
  setwd(direct)
  files.list <- list.files()
  files.total <- length(files.list)
  
  # Variable to keep track of the count that matches the conditions
  modifiedCount <- 1 
  
  
  # Open file for reading in binary mode
  raw.file <- file(files.list[1], "rb")
  
  # Read file header
  raw.header <- readBin(raw.file, integer(), 1)
  raw.totalpulses <- readBin(raw.file, integer(), 1)
  raw.blank <- readBin(raw.file, integer(), 1)
  raw.blank2 <- readBin(raw.file, integer(), 1)
  raw.acqTime <- readBin(raw.file, integer(), 1)
  raw.headerRemaining.bin <- readBin(raw.file, raw(), 4100)
  
  # Create empty vectors for storing pulse information
  # This will work much faster than expanding them on the fly.
  c.beginTimes <- rep(0, raw.totalpulses)
  c.heights <- rep(0, raw.totalpulses)
  c.coincs <- rep(0, raw.totalpulses)
  c.widths <- rep(0, raw.totalpulses)
  q
  
  
  for (j in seq(1, raw.totalpulses)) {
    
    # print(j/raw.totalpulses)
    # We don't need a progress bar here. It's very fast.
    
    pulse.beginTime <- readBin(raw.file, integer(), n=1, size=8)
    pulse.width <- readBin(raw.file, integer(), n=1, size=2)
    pulse.coinc <- readBin(raw.file, integer(), n=1, size=2)
    pulse.height <- readBin(raw.file, integer(), n=1, size=4)
    
    if (pulse.width < WidthMax && pulse.width > WidthMin && 
        pulse.height > HeightMin&& pulse.height<HeightMax){
      c.beginTimes[modifiedCount] <- pulse.beginTime
      c.heights[modifiedCount] <- pulse.height
      c.coincs[modifiedCount] <- pulse.coinc
      c.widths[modifiedCount] <- pulse.width
      modifiedCount <- modifiedCount + 1
      
      
      
    }
    
  }
  
  # Trim the array
  c.beginTimes <- c.beginTimes[1:modifiedCount-1]
  c.heights <- c.heights[1:modifiedCount-1]
  c.coincs <- c.coincs[1:modifiedCount-1]
  c.widths <- c.widths[1:modifiedCount-1]
  
  
  # Organize this data into a data frame
  pulsesC.data <- data.frame(
    pulseId = seq(1, modifiedCount-1),
    pulseBeginTime = c.beginTimes,
    pulseHeight = c.heights,
    pulseCoincidence = c.coincs,
    pulseWidth = c.widths
  )
  
  # Normalize the raw numbers
  pulsesC.data$pulseHeightNorm <- pulsesC.data$pulseHeight*100/2147483647
  pulsesC.data$pulseWidthMs <- pulsesC.data$pulseWidth/48
  
  # Produce plot
  # Heat map of Width vs Height and counts
  #   dir <- paste(direct,'1',sep="")
  #   dir <- paste(dir,"M.jpeg",sep="")
  #   jpeg(file=dir)
  #   #   
  #   p <- ggplot(pulsesC.data, aes(x=pulseHeightNorm, y=pulseWidthMs))
  #   p + geom_point()
  #   #   
  #   print(
  #     
  #     p + theme_bw() + stat_bin2d(drop=F, binwidth=c(0.04, 0.06)) + 
  #       scale_x_continuous("pulse height (a.u.)", expand=c(0,0), limits=c(1.5,10)) +
  #       scale_y_continuous("pulse width (ms)", expand=c(0,0), limits=c(0,15)) + 
  #       coord_fixed(ratio=0.5) +
  #       scale_fill_gradient(
  #         trans="log10", na.value="#000000", 
  #         low="#000000", high="#71ff71",
  #         limits=c(1,10)
  #       ) +
  #       theme(
  #         panel.grid = element_blank()   
  #       )
  #   )
  #dev.off()
  # Returns control data
  return(pulsesC.data)
}




# Maximum and minimum boundarys for the signals to be processed, and factor to enhance control subtraction 
# height.maxthresh <- .Machine$integer.max
# height.minthresh <- 0

height.maxthresh <- as.numeric(readline("Enter upper threshold for heights: "))*2147483647/100
height.minthresh <- as.numeric(readline("Enter lower threshold for heights: "))*2147483647/100
width.maxthresh <- as.numeric(readline("Enter upper threshold for widths: "))*48
width.minthresh <- as.numeric(readline("Enter lower threshold for widths: "))*48
alpha <- as.numeric(readline("Enter 1 for Individual plots, 2 for no plots: "))
beta <- as.numeric(readline("Enter 3 for csv files, 4 for no csv: "))

# Calls function that runs on control file and return necessary data
Control.data <- calculateControl(ControlDirectory,height.maxthresh,height.minthresh,width.maxthresh,width.minthresh)

# Calls function that produces 3 plots of test data with control data subtracted, also returns vector of integrals, each index represents
# one file processe
Result<- calculateandplotPulseData(TestDirectory,Control.data,height.maxthresh,height.minthresh,width.maxthresh,width.minthresh,alpha,beta)

# Need to get standard curve to find line of best fit
# Plots the cell count vs time with error bars
