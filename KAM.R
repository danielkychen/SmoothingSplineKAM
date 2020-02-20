#Calculation of KAM while changing step width and FPA

rm(list=ls())

library(pracma)
library(signal)
library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(rgl)
library(car)
library(PerformanceAnalytics)
library(earth)
library(plot3D)

#Functions for script-------------------------------------------------------------------------------------------------------------------

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Calculates midpoint of two markers
midpoint <- function(marker1, marker2){
  midpoint <- mean(c(marker1,marker2))
  
  return(midpoint)
}

#Position Vector calculation based on two markers
positionVec <- function(origin, endPoint){
  
  posVec <- endPoint - origin
  
  return(posVec)
}

#Knee Adduction Moment calculation
KAM <- function(posVec, GRF){
  KAM <- cross(posVec, GRF)
  
  return(KAM)
}

#Filter signal
filterSignal <- function(signal, cutoff, type, sampleFreq){
  
  bf <- butter(4, (2*cutoff)/(sampleFreq), type = type)
  y <- filtfilt(bf, signal)
  
  return(y)
}

#Calculates COP from force plate data
#(238.5, 812.5, 0) for force plate 1, (727.5, 812.5, 0) for force plate 2??
COP <- function(dataframe, leftorright = "left"){
  attach(dataframe)
  h <- 15 #height of padding
  
  if (leftorright == "left"){
    COPx <- ((-h*Force.Plate.1_Force.X_N - Force.Plate.1_Moment.Y_Nmm)/Force.Plate.1_Force.Z_N) + Force.Plate.1_Ref.X_mm
    COPy <- ((h*Force.Plate.1_Force.Y_N + Force.Plate.1_Moment.X_Nmm)/Force.Plate.1_Force.Z_N) + Force.Plate.1_Ref.Y_mm
    
  } else if (leftorright == "right"){
    COPx <- ((-h*Force.Plate.2_Force.X_N - Force.Plate.2_Moment.Y_Nmm)/Force.Plate.2_Force.Z_N) + Force.Plate.2_Ref.X_mm
    COPy <- ((h*Force.Plate.2_Force.Y_N + Force.Plate.2_Moment.X_Nmm)/Force.Plate.2_Force.Z_N) + Force.Plate.2_Ref.Y_mm
  }
  
  
  # COPx <- ((-orgZ*Fx - My)/Fz) + orgX
  # COPy <- ((-orgZ*Fy + Mx)/Fz) + orgY
  COPz <- 0
  
  # if (component == "x"){
  #   COP = COPx
  # } else if (component == "y"){
  #   COP = COPy
  # } else if (component == "z"){
  #   COP = COPz
  # } else {
  #   COP = "Error"
  # }
  detach(dataframe)
  COP <- list("x" = COPx, "y" = COPy, "z" = COPz)
  
  return(COP)
}

#Calculate magnitude of vector
mag <- function(x, y, z){
  mag <- sqrt((x*x)+(y*y)+(z*z))
  return(mag)
}

#Rotate vector about Global Z Axis
rotate <- function(theta, x, y, z){
  x2 <- x*cos(theta*pi/180) + y*sin(theta*pi/180)
  y2 <- -x*sin(theta*pi/180) + y*cos(theta*pi/180)
  
  return(c(x2, y2, z))
}

projection <- function(x, y, z, x2, y2, z2){
  #determine unit vector of x2, y2, z2
  Y <- c(x2, y2, z2)/(mag(x2, y2, z2))
  M <- c(x, y, z)
  
  MdotY <- dot(M, Y)
  
  Y2 <- MdotY*(Y)
  
  return(Y2)
}

filterDataframe <- function(dataframe){
  
  #Replace NAs with 0
  dataframe[is.na(dataframe)] <- 0
  
  for (i in 1:ncol(dataframe)){
    dataframe[,i] <- with(data = dataframe, expr = filterSignal(dataframe[,i], cutoff = 10, type = 'low', sampleFreq = 100))
  }
  
  return(dataframe)
}

extractVariables <- function(dataframe, variables){
  #variablesToKeep <- colnames(dataframe)
  keep <- rep(FALSE, length(colnames(dataframe)))
  
  for (i in 1:length(variables)){
    for (j in 1:length(colnames(dataframe))){
      
      if ((strsplit(colnames(dataframe)[j], split = "_")[[1]][1]) == variables[i]){
        keep[j] <- TRUE
      }
    }
  }
  
  return(dataframe[,keep])
} 

#Calculate and update dataframe with FPA
FPA <- function(dataframe){
  LeftDot <- c()
  RightDot <- c()
  
  #Calculate Calcaneal to Toe vector
  LeftFootVec_X <- positionVec(origin = dataframe$LC_X_mm, endPoint = dataframe$LT_X_mm)
  RightFootVec_X <- positionVec(origin = dataframe$RC_X_mm, endPoint = dataframe$RT_X_mm)
  
  LeftFootVec_Y <- positionVec(origin = dataframe$LC_Y_mm, endPoint = dataframe$LT_Y_mm)
  RightFootVec_Y <- positionVec(origin = dataframe$RC_Y_mm, endPoint = dataframe$RT_Y_mm)
  
  #Y Axis vector with which angle is calculated with
  yAxis_X <- 0
  yAxis_Y <- 1
  
  yAxis <- c(0, 1)
  
  #Calculate Angle between vectors
  for (i in 1:nrow(dataframe)){
    LeftDot[i] <- dot(c(LeftFootVec_X[i], LeftFootVec_Y[i]), yAxis)
    RightDot[i] <- dot(c(RightFootVec_X[i], RightFootVec_Y[i]), yAxis)
  }
  
  LeftMag <- mag(x = LeftFootVec_X, y = LeftFootVec_Y, z = 0)
  RightMag <- mag(x = RightFootVec_X, y = RightFootVec_Y, z = 0)
  
  LeftAngle <- 180*(acos(LeftDot/(LeftMag)))/pi
  RightAngle <- 180*(acos(RightDot/(RightMag)))/pi
  
  #Correct for Toe-in and Toe-out (-ve and +ve)
  for (i in 1:length(LeftAngle)){
    if(LeftFootVec_X[i] > 0){
      LeftAngle[i] <- LeftAngle[i]*(-1)
    }
    
    if(RightFootVec_X[i] < 0){
      RightAngle[i] <- RightAngle[i]*(-1)
    }
    
  }
    
  Angle <- list("Left" = LeftAngle, "Right" = RightAngle)
    
  return(Angle)
}

StepWidth <- function(dataframe){
  
  LeftHeel <- rep(0, times = nrow(dataframe))
  RightHeel <- rep(0, times = nrow(dataframe)) #dataframe$RC_X_mm
  width <- rep(0, times = nrow(dataframe))
  
  for (i in 1:nrow(dataframe)){
    if (dataframe[i,'StartLeftStep']){
      LeftHeel[i] <- dataframe[i,'LC_X_mm']
    } else {
        if (i != 1){
          LeftHeel[i] <- LeftHeel[i-1]
        }
    }
    
    if (dataframe[i,'StartRightStep']){
      RightHeel[i] <- dataframe[i,'RC_X_mm']
    } else {
        if (i != 1){
          RightHeel[i] <- RightHeel[i-1]
        }
    }
    
    width[i] <- RightHeel[i] - LeftHeel[i]
  }
  
  return(width)
}

Step <- function(dataframe){
    LeftStep <- rep(FALSE, times = nrow(dataframe))
    RightStep <- rep(FALSE, times = nrow(dataframe))
    StartLeftStep <- rep(FALSE, times = nrow(dataframe))
    StartRightStep <- rep(FALSE, times = nrow(dataframe))
    EndLeftStep <- rep(FALSE, times = nrow(dataframe))
    EndRightStep <- rep(FALSE, times = nrow(dataframe))
    
  for (i in 1:nrow(dataframe)){
    if ((dataframe[i,'Force.Plate.1_Force.Z_N'] > 80)){
      LeftStep[i] <- TRUE
    }
    
    if ((dataframe[i,'Force.Plate.2_Force.Z_N'] > 80)){
      RightStep[i] <- TRUE
    }
    
    if ((i != 1) && (LeftStep[i] > LeftStep[i-1])){
      StartLeftStep[i] <- TRUE
    }
    
    if ((i != 1) && (RightStep[i] > RightStep[i-1])){
      StartRightStep[i] <- TRUE
    }
    
    if ((i != 1) && (LeftStep[i] < LeftStep[i-1])){
      EndLeftStep[i] <- TRUE
    }
    
    if ((i != 1) && (RightStep[i] < RightStep[i-1])){
      EndRightStep[i] <- TRUE
    }
  }
    
  stepList <- list("StartLeft" = StartLeftStep, "StartRight" = StartRightStep, "EndLeft" = EndLeftStep, "EndRight" = EndRightStep)
  
  return(stepList)
}

CalculateKam <- function(dataframe){
  
  headers <- colnames(dataframe)
  
  for (i in 1:length(headers)){
    if ((headers[i] == "LTibiaToGlobal_Z_mm")||(headers[i] == "LTibiaToGlobal_Z_degrees")||(headers[i] == "LTibiaToGlobal_Z")){
      headers[i] <- "LTibiaToGlobal_Z_deg"
    }
    
    if ((headers[i] == "RTibiaToGlobal_Z_mm")||(headers[i] == "RTibiaToGlobal_Z_degrees")||(headers[i] == "RTibiaToGlobal_Z")){
      headers[i] <- "RTibiaToGlobal_Z_deg"
    }
  }
  colnames(dataframe) <- headers
  
  
  
  attach(dataframe)
  
  
  KAMx_Left <- c()
  KAMy_Left <- c()
  KAMz_Left <- c()
  KAMx_Right <- c()
  KAMy_Right <- c()
  KAMz_Right <- c()
  LTibKAMx <- c()
  LTibKAMy <- c()
  LTibKAMz <- c()
  RTibKAMx <- c()
  RTibKAMy <- c()
  RTibKAMz <- c()
  
  LKAM <- c()
  RKAM <- c()
  #UNITS ARE Nmm
  
  for (i in 1:nrow(dataframe)){
    if (Force.Plate.1_Force.Z_N[i] > 80){
      KAMx_Left <- append(KAMx_Left, KAM(c(LposX[i], LposY[i], LposZ[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[1])
      KAMy_Left <- append(KAMy_Left, KAM(c(LposX[i], LposY[i], LposZ[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[2])
      KAMz_Left <- append(KAMz_Left, KAM(c(LposX[i], LposY[i], LposZ[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[3])
      
      #Rotate KAM to Tibia Frame
      LTibKAMx <- append(LTibKAMx, rotate(LTibiaToGlobal_Z_deg[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[1])
      LTibKAMy <- append(LTibKAMy, rotate(LTibiaToGlobal_Z_deg[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[2])
      LTibKAMz <- append(LTibKAMz, rotate(LTibiaToGlobal_Z_deg[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[3])
      
    } else {
      LTibKAMx <- append(LTibKAMx, 0)
      LTibKAMy <- append(LTibKAMy, 0)
      LTibKAMz <- append(LTibKAMz, 0)
      
      KAMx_Left <- append(KAMx_Left, 0)
      KAMy_Left <- append(KAMy_Left, 0)
      KAMz_Left <- append(KAMz_Left, 0)
    }
    
    if (Force.Plate.2_Force.Z_N[i] > 80){
      KAMx_Right <- append(KAMx_Right, KAM(c(RposX[i], RposY[i], RposZ[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[1])
      KAMy_Right <- append(KAMy_Right, KAM(c(RposX[i], RposY[i], RposZ[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[2])
      KAMz_Right <- append(KAMz_Right, KAM(c(RposX[i], RposY[i], RposZ[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[3])
      
      #Rotate KAM to Tibia Frame
      RTibKAMx <- append(RTibKAMx, rotate(RTibiaToGlobal_Z_deg[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[1])
      RTibKAMy <- append(RTibKAMy, rotate(RTibiaToGlobal_Z_deg[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[2])
      RTibKAMz <- append(RTibKAMz, rotate(RTibiaToGlobal_Z_deg[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[3])
      
    } else {
      KAMx_Right <- append(KAMx_Right, 0)
      KAMy_Right <- append(KAMy_Right, 0)
      KAMz_Right <- append(KAMz_Right, 0)
      
      RTibKAMx <- append(RTibKAMx, 0)
      RTibKAMy <- append(RTibKAMy, 0)
      RTibKAMz <- append(RTibKAMz, 0)
    }
  }
  
  detach(dataframe)
  
  LKAM <- mag(LTibKAMx, LTibKAMy, LTibKAMz)
  RKAM <- mag(RTibKAMx, RTibKAMy, RTibKAMz)
  
  TibList <- list("Left" = LKAM, "Right" = RKAM)
  
  return(TibList)
}

PosVec <- function(dataframe){
  
  LposX <- mapply(positionVec, dataframe$LKJC_X_mm, dataframe$Force.Plate.1_COP.X_mm)
  LposY <- mapply(positionVec, dataframe$LKJC_Y_mm, dataframe$Force.Plate.1_COP.Y_mm)
  LposZ <- mapply(positionVec, dataframe$LKJC_Z_mm, dataframe$Force.Plate.1_COP.Z_mm)
  
  RposX <- mapply(positionVec, dataframe$RKJC_X_mm, dataframe$Force.Plate.2_COP.X_mm)
  RposY <- mapply(positionVec, dataframe$RKJC_Y_mm, dataframe$Force.Plate.2_COP.Y_mm)
  RposZ <- mapply(positionVec, dataframe$RKJC_Z_mm, dataframe$Force.Plate.2_COP.Z_mm)
  
  PosList <- list("Lx" = LposX, "Ly" = LposY, "Lz" = LposZ, "Rx" = RposX, "Ry" = RposY, "Rz" = RposZ)
  
  return(PosList)
}

ScatterPlotInterative <- function(dataframe, leftorright = "left"){
  attach(dataframe)
  if (leftorright == "left"){
    plot <- plot3d(FPA_Left, SW, LKAM_norm, col="red", size=3)
  } else {
    plot <- plot3d(FPA_Right, SW, RKAM_norm, col="red", size=3)
  }
  detach(dataframe)
  
  return(plot)
}

AllScatter <- function(listofDataFrames, type = "multi"){
  AllData <- data.frame()
  for (i in 1:length(listofDataFrames)){
    AllData <- rbind.data.frame(AllData, listofDataFrames[[i]])
  }
  if (type == "multi"){
    plot <- chart.Correlation(AllData, histogram = FALSE)
  } else if (type == "3d"){
    plot <- plot3d(AllData$secondPeakFPA, AllData$secondPeakSW, AllData$maxPeakKAM, col="red", size=3)
  }
  
  
  return(plot)
}

segmentStep <- function(dataframe, leftorright = "left"){
  attach(dataframe)
  
  count <- 0
  startIndex <- 0
  endIndex <- 0
  midIndex <- 0
  FPA <- c()
  sw <- c()
  peakKAM <- c()
  firstKAM <- c()
  secondKAM <- c()

  
  stepKAMarray <- c()
  
  if (leftorright == "left"){
    for (i in 1:nrow(dataframe)){
      if (StartLeftStep[i]){
        startIndex <- i
      }
      
      if (EndLeftStep[i]){
        endIndex <- i
        midIndex <- ((endIndex - startIndex)/2) + startIndex
        FPA <- append(FPA, FPA_Left[midIndex])
        sw <- append(sw, SW[midIndex])
        stepKAMarray <- c(startIndex:endIndex)
        peakKAM <- append(peakKAM, max(LKAM_norm[stepKAMarray]))
      }
    }
  } else if (leftorright == "right"){
    for (i in 1:nrow(dataframe)){
      if (StartRightStep[i]){
        startIndex <- i
      }
      
      if (EndRightStep[i]){
        endIndex <- i
        midIndex <- ((endIndex - startIndex)/2) + startIndex
        FPA <- append(FPA, FPA_Left[midIndex])
        sw <- append(sw, SW[midIndex])
        stepKAMarray <- c(startIndex:endIndex)
        peakKAM <- append(peakKAM, max(LKAM_norm[stepKAMarray]))
      }
    }
  }
  
  #Remove first and second steps because not legit
  new_dataframe <- data.frame(FPA[c(-1,-2,-3)], sw[c(-1,-2,-3)], peakKAM[c(-1,-2,-3)])
  colnames(new_dataframe) <- c("FPA", "Step Width", "Peak KAM")
  
  new_dataframe <- new_dataframe[new_dataframe$`Peak KAM`!=0,]
  
  detach(dataframe)
  
  return(new_dataframe)
}

plotSteps <- function(dataframe, leftorright = "left"){
  dataframe <- cbind.data.frame(dataframe, "Time" = 1:nrow(dataframe))
  
  if (leftorright == "left"){
    Step_dataframe <- dataframe[,c("Time", "StartLeftStep", "EndLeftStep", "LKAM_norm")]
  } else if (leftorright == "right"){
    Step_dataframe <- dataframe[,c("Time", "StartRightStep", "EndRightStep", "RKAM_norm")]
  }
  
  melt_dataframe <- melt(data = Step_dataframe, id.vars = "Time") #id.vars is usually the x axis
  Plot <- ggplot(data = melt_dataframe, aes(x = melt_dataframe$Time, y = melt_dataframe$value, colour = melt_dataframe$variable)) + geom_line()
  
  return(Plot)
}

clipAt <- function(dataframe, leftorright = "left", clipStart = 1, clipEnd = 1){
  attach(dataframe)
  StartCount <- 0
  EndCount <- 0
  
  if (leftorright == "left"){
    totalEnd <- sum(dataframe$EndLeftStep)
    for (i in 1:nrow(dataframe)){
      if (dataframe[i, "StartLeftStep"]){
        StartCount <- StartCount + 1
        
        if (StartCount == clipStart){
          dataframe[(1:i), "LKAM_norm"] <- 0
          dataframe[(1:(i-1)), c("StartLeftStep","EndLeftStep")] <- FALSE
        }
      }
      
      if (dataframe[i, "EndLeftStep"]){
        EndCount <- EndCount + 1
        
        if (EndCount == totalEnd - (clipEnd - 1)){
          dataframe[(i:nrow(dataframe)), "LKAM_norm"] <- 0
          dataframe[((i+1):nrow(dataframe)), c("StartLeftStep","EndLeftStep")] <- FALSE
          break
        }
      }
    }
  } else if (leftorright == "right"){
    totalEnd <- sum(dataframe$EndRightStep)
    for (i in 1:nrow(dataframe)){
      if (dataframe[i, "StartRightStep"]){
        StartCount <- StartCount + 1
        
        if (StartCount == clipStart){
          dataframe[(1:i), "RKAM_norm"] <- 0
          dataframe[(1:(i-1)), c("StartRightStep","EndRightStep")] <- FALSE
        }
      }
      
      if (dataframe[i, "EndRightStep"]){
        EndCount <- EndCount + 1
        
        if (EndCount == totalEnd - (clipEnd - 1)){
          dataframe[(i:nrow(dataframe)), "RKAM_norm"] <- 0
          dataframe[((i+1):nrow(dataframe)), c("StartRightStep","EndRightStep")] <- FALSE
          break
        }
      }
    }
  }
  detach(dataframe)
  return(dataframe)
}

NearestPeak <- function(dataframe, peaks, leftorright = "left", number = 3){
  firstPeak <- peaks[1]
  secondPeak <- peaks[2]
  
  firstRange <- c((firstPeak - number):(firstPeak + number))
  secondRange <- c((secondPeak - number):(secondPeak + number))
  
  if (leftorright == "left"){
    for (i in 1:length(firstRange)){
      if (max(dataframe[firstRange, "LKAM_norm"]) == (dataframe[firstRange[i], "LKAM_norm"])){
        firstPeak <- firstRange[i]
        break
      }
    }
    
    for (i in 1:length(secondRange)){
      if (max(dataframe[secondRange, "LKAM_norm"]) == (dataframe[secondRange[i], "LKAM_norm"])){
        secondPeak <- secondRange[i]
        break
      }
    }
  } else if (leftorright == "right"){
    for (i in 1:length(firstRange)){
      if (max(dataframe[firstRange, "RKAM_norm"]) == (dataframe[firstRange[i], "RKAM_norm"])){
        firstPeak <- firstRange[i]
        break
      }
    }
    
    for (i in 1:length(secondRange)){
      if (max(dataframe[secondRange, "RKAM_norm"]) == (dataframe[secondRange[i], "RKAM_norm"])){
        secondPeak <- secondRange[i]
        break
      }
    }
  }
  
  
  peaker <- c(firstPeak, secondPeak)
  return(peaker)
}

PeakSelect <- function(dataframe, leftorright = "left"){
  attach(dataframe)
  startIndex <- 0
  endIndex <- 0
  peaks <- c()
  firstPeak <- c()
  secondPeak <- c()
  par(mfrow=c(1,1))
  
  if (leftorright == "left"){
    for (i in 1:nrow(dataframe)){
      if (StartLeftStep[i]){
        startIndex <- i
      }
      
      if (EndLeftStep[i]){
        endIndex <- i
        print(plot(x = startIndex:endIndex, y = LKAM_norm[startIndex:endIndex], type = 'l'))
        # print(plot(x = startIndex:endIndex, y = Force.Plate.1_COP.X_mm[startIndex:endIndex], type = 'l'))
        # print(plot(x = startIndex:endIndex, y = Force.Plate.1_COP.Y_mm[startIndex:endIndex], type = 'l'))
        # print(plot(x = startIndex:endIndex, y = Force.Plate.1_Force.Z_N[startIndex:endIndex], type = 'l'))
        peaks <- identify(x = startIndex:endIndex, y = LKAM_norm[startIndex:endIndex], n = 2)
        if (length(peaks) == 2){
          peaker <- NearestPeak(dataframe, peaks+startIndex-1, "left")
        } else {peaker <- peaks}
        firstPeak <- append(firstPeak, peaker[1])
        secondPeak <- append(secondPeak, peaker[2])
      }
    }
  } else if (leftorright == "right"){
    for (i in 1:nrow(dataframe)){
      if (StartRightStep[i]){
        startIndex <- i
      }
      
      if (EndRightStep[i]){
        endIndex <- i
        #peakFinder here
        print(plot(x = startIndex:endIndex, y = RKAM_norm[startIndex:endIndex], type = 'l'))
        # print(plot(x = startIndex:endIndex, y = Force.Plate.2_COP.X_mm[startIndex:endIndex], type = 'l'))
        # print(plot(x = startIndex:endIndex, y = Force.Plate.2_COP.Y_mm[startIndex:endIndex], type = 'l'))
        # print(plot(x = startIndex:endIndex, y = Force.Plate.2_Force.Z_N[startIndex:endIndex], type = 'l'))
        peaks <- identify(x = startIndex:endIndex, y = RKAM_norm[startIndex:endIndex], n = 2)
        if (length(peaks) == 2){
          peaker <- NearestPeak(dataframe, peaks+startIndex-1, "right")
        } else {peaker <- peaks}
        firstPeak <- append(firstPeak, peaker[1])
        secondPeak <- append(secondPeak, peaker[2])
      }
    }
  }
  
  peakers <- list("first" = firstPeak, "second" = secondPeak)
  
  return(peakers)
}

makeZero <- function(dataframe){
  attach(dataframe)
  
    for (i in 1:nrow(dataframe)){
      if (Force.Plate.1_Force.Z_N[i] < 80){
        Force.Plate.1_Force.X_N[i] <- 0
        Force.Plate.1_Force.Y_N[i] <- 0
        Force.Plate.1_Force.Z_N[i] <- 0
        Force.Plate.1_Moment.X_Nmm[i] <- 0
        Force.Plate.1_Moment.Y_Nmm[i] <- 0
        Force.Plate.1_Moment.Z_Nmm[i] <- 0
      }
      if (Force.Plate.2_Force.Z_N[i] < 80){
        Force.Plate.2_Force.X_N[i] <- 0
        Force.Plate.2_Force.Y_N[i] <- 0
        Force.Plate.2_Force.Z_N[i] <- 0
        Force.Plate.2_Moment.X_Nmm[i] <- 0
        Force.Plate.2_Moment.Y_Nmm[i] <- 0
        Force.Plate.2_Moment.Z_Nmm[i] <- 0
      }
    }
    dataframe$Force.Plate.1_Force.X_N <- Force.Plate.1_Force.X_N
    dataframe$Force.Plate.1_Force.Y_N <- Force.Plate.1_Force.Y_N
    dataframe$Force.Plate.1_Force.Z_N <- Force.Plate.1_Force.Z_N
    dataframe$Force.Plate.1_Moment.X_Nmm <- Force.Plate.1_Moment.X_Nmm
    dataframe$Force.Plate.1_Moment.Y_Nmm <- Force.Plate.1_Moment.Y_Nmm
    dataframe$Force.Plate.1_Moment.Z_Nmm <- Force.Plate.1_Moment.Z_Nmm

    dataframe$Force.Plate.2_Force.X_N <- Force.Plate.2_Force.X_N
    dataframe$Force.Plate.2_Force.Y_N <- Force.Plate.2_Force.Y_N
    dataframe$Force.Plate.2_Force.Z_N <- Force.Plate.2_Force.Z_N
    dataframe$Force.Plate.2_Moment.X_Nmm <- Force.Plate.2_Moment.X_Nmm
    dataframe$Force.Plate.2_Moment.Y_Nmm <- Force.Plate.2_Moment.Y_Nmm
    dataframe$Force.Plate.2_Moment.Z_Nmm <- Force.Plate.2_Moment.Z_Nmm
  
  detach(dataframe)
  
  return(dataframe)
}

peakDataframe <- function(dataframe, peakvalues, leftorright = "left"){

  firstPeak <- peakvalues$first
  secondPeak <- peakvalues$second
  firstPeak <- firstPeak[!is.na(firstPeak)]
  secondPeak <- secondPeak[!is.na(secondPeak)]
  
  if (leftorright == "left"){
    FPA <- dataframe[firstPeak, "FPA_Left"]
    firstPeakKAM <- dataframe[firstPeak, "LKAM_norm"]
    secondPeakKAM <- dataframe[secondPeak, "LKAM_norm"]
  } else if (leftorright == "right"){
    FPA <- dataframe[firstPeak, "FPA_Right"]
    firstPeakKAM <- dataframe[firstPeak, "RKAM_norm"]
    secondPeakKAM <- dataframe[secondPeak, "RKAM_norm"]
  }
  
  StepWidth <- dataframe[firstPeak, "SW"]
  maxPeakKAM <- c()
  
  for (i in 1:length(firstPeak)){
    if (firstPeakKAM[i] > secondPeakKAM[i]){
      maxPeakKAM[i] <- firstPeakKAM[i]
    } else {
      maxPeakKAM[i] <- secondPeakKAM[i]
    }
  }
  
  new_dataframe <- data.frame(FPA, StepWidth, firstPeakKAM,
                              secondPeakKAM, maxPeakKAM)
  
  return(new_dataframe)
}

#---------------------------------------------------------------------------------------------------------------------------------------

dir <- setwd(dir = 'C:\\Users\\bioengsu\\Desktop\\KAM\\Cleaned CSV\\Trial 1')

subjects <- list.files(dir)

#Read data from subjects----------------------------------------------------------------------------------------------------------------
AndrewBollen.df <- read.csv(file = subjects[1], header = TRUE, stringsAsFactors = FALSE)
AndrewWong.df <- read.csv(file = subjects[2], header = TRUE, stringsAsFactors = FALSE)
Daniel.df <- read.csv(file = subjects[3], header = TRUE, stringsAsFactors = FALSE)
Frank.df <- read.csv(file = subjects[4], header = TRUE, stringsAsFactors = FALSE)
Jasmine.df <- read.csv(file = subjects[5], header = TRUE, stringsAsFactors = FALSE)
Jonathan.df <- read.csv(file = subjects[6], header = TRUE, stringsAsFactors = FALSE)
Kelly.df <- read.csv(file = subjects[7], header = TRUE, stringsAsFactors = FALSE)
Phoebe.df <- read.csv(file = subjects[8], header = TRUE, stringsAsFactors = FALSE)
Rakesh.df <- read.csv(file = subjects[9], header = TRUE, stringsAsFactors = FALSE)
Thor.df <- read.csv(file = subjects[10], header = TRUE, stringsAsFactors = FALSE)
Wayne.df <- read.csv(file = subjects[11], header = TRUE, stringsAsFactors = FALSE)
#-------------------------------------------------------------------------------------------------------------

#Keep only Variables of interest------------------------------------------------------------------------------
#RC, RT, LC, LT, RTIBTOGLOBAL_Z, LTIBTOGLOBAL_Z, RKJC, LKJC, FORCE, MOMENT
variables <- c("RC", "RT", "LC", "LT", "RKJC", "LKJC", "RTibiaToGlobal", "LTibiaToGlobal", "Force.Plate.1", "Force.Plate.2")
  
AndrewBollen.df <- extractVariables(AndrewBollen.df, variables = variables)
AndrewWong.df <- extractVariables(AndrewWong.df, variables = variables)
Daniel.df <- extractVariables(Daniel.df, variables = variables)
Frank.df <- extractVariables(Frank.df, variables = variables)
Jasmine.df <- extractVariables(Jasmine.df, variables = variables)
Jonathan.df <- extractVariables(Jonathan.df, variables = variables)
Kelly.df <- extractVariables(Kelly.df, variables = variables)
Phoebe.df <- extractVariables(Phoebe.df, variables = variables)
Rakesh.df <- extractVariables(Rakesh.df, variables = variables)
Thor.df <- extractVariables(Thor.df, variables = variables)
Wayne.df <- extractVariables(Wayne.df, variables = variables)
#-------------------------------------------------------------------------------------------------------------

#Filter Force, moment, marker data (if not done in Vicon)------------------------------------------------------------------------------------------------------------
# AndrewBollen.df <- filterDataframe(AndrewBollen.df)
# AndrewWong.df <- filterDataframe(AndrewWong.df)
# Daniel.df <- filterDataframe(Daniel.df)
# Frank.df <- filterDataframe(Frank.df)
# Jasmine.df <- filterDataframe(Jasmine.df)
# Jonathan.df <- filterDataframe(Jonathan.df)
# Kelly.df <- filterDataframe(Kelly.df)
# Phoebe.df <- filterDataframe(Phoebe.df)
# Rakesh.df <- filterDataframe(Rakesh.df)
# Thor.df <- filterDataframe(Thor.df)
# Wayne.df <- filterDataframe(Wayne.df)
#----------------------------------------------------------------------------------------------------------------------------------------

#Multiple Force plate Y axis by -1 ------------------------------------
AndrewBollen.df$Force.Plate.1_Force.Y_N <- AndrewBollen.df$Force.Plate.1_Force.Y_N * (-1)
AndrewWong.df$Force.Plate.1_Force.Y_N <- AndrewWong.df$Force.Plate.1_Force.Y_N * (-1)
Daniel.df$Force.Plate.1_Force.Y_N <- Daniel.df$Force.Plate.1_Force.Y_N * (-1)
Frank.df$Force.Plate.1_Force.Y_N <- Frank.df$Force.Plate.1_Force.Y_N * (-1)
Jasmine.df$Force.Plate.1_Force.Y_N <- Jasmine.df$Force.Plate.1_Force.Y_N * (-1)
Jonathan.df$Force.Plate.1_Force.Y_N <- Jonathan.df$Force.Plate.1_Force.Y_N * (-1)
Kelly.df$Force.Plate.1_Force.Y_N <- Kelly.df$Force.Plate.1_Force.Y_N * (-1)
Phoebe.df$Force.Plate.1_Force.Y_N <- Phoebe.df$Force.Plate.1_Force.Y_N * (-1)
Rakesh.df$Force.Plate.1_Force.Y_N <- Rakesh.df$Force.Plate.1_Force.Y_N * (-1)
Thor.df$Force.Plate.1_Force.Y_N <- Thor.df$Force.Plate.1_Force.Y_N * (-1)
Wayne.df$Force.Plate.1_Force.Y_N <- Wayne.df$Force.Plate.1_Force.Y_N * (-1)

AndrewBollen.df$Force.Plate.2_Force.Y_N <- AndrewBollen.df$Force.Plate.2_Force.Y_N * (-1)
AndrewWong.df$Force.Plate.2_Force.Y_N <- AndrewWong.df$Force.Plate.2_Force.Y_N * (-1)
Daniel.df$Force.Plate.2_Force.Y_N <- Daniel.df$Force.Plate.2_Force.Y_N * (-1)
Frank.df$Force.Plate.2_Force.Y_N <- Frank.df$Force.Plate.2_Force.Y_N * (-1)
Jasmine.df$Force.Plate.2_Force.Y_N <- Jasmine.df$Force.Plate.2_Force.Y_N * (-1)
Jonathan.df$Force.Plate.2_Force.Y_N <- Jonathan.df$Force.Plate.2_Force.Y_N * (-1)
Kelly.df$Force.Plate.2_Force.Y_N <- Kelly.df$Force.Plate.2_Force.Y_N * (-1)
Phoebe.df$Force.Plate.2_Force.Y_N <- Phoebe.df$Force.Plate.2_Force.Y_N * (-1)
Rakesh.df$Force.Plate.2_Force.Y_N <- Rakesh.df$Force.Plate.2_Force.Y_N * (-1)
Thor.df$Force.Plate.2_Force.Y_N <- Thor.df$Force.Plate.2_Force.Y_N * (-1)
Wayne.df$Force.Plate.2_Force.Y_N <- Wayne.df$Force.Plate.2_Force.Y_N * (-1)
#----------------------------------------------------------------------

#Make Zero if under 80N threshold for Fz-------------------------------
AndrewBollen.df <- makeZero(AndrewBollen.df)
AndrewWong.df <- makeZero(AndrewWong.df)
Daniel.df <- makeZero(Daniel.df)
Frank.df <- makeZero(Frank.df)
Jasmine.df <- makeZero(Jasmine.df)
Jonathan.df <- makeZero(Jonathan.df)
Kelly.df <- makeZero(Kelly.df)
Phoebe.df <- makeZero(Phoebe.df)
Rakesh.df <- makeZero(Rakesh.df)
Thor.df <- makeZero(Thor.df)
Wayne.df <- makeZero(Wayne.df)
#----------------------------------------------------------------------

#COP Calculation------------------------------------------------------------------------------------
#Overwrite the COP values from Vicon as they were based off pre-filtered data
AndrewBollen.df$Force.Plate.1_COP.X_mm <- COP(AndrewBollen.df, leftorright = "left")$x
AndrewBollen.df$Force.Plate.1_COP.Y_mm <- COP(AndrewBollen.df, leftorright = "left")$y

AndrewBollen.df$Force.Plate.2_COP.X_mm <- COP(AndrewBollen.df, leftorright = "right")$x
AndrewBollen.df$Force.Plate.2_COP.Y_mm <- COP(AndrewBollen.df, leftorright = "right")$y

AndrewWong.df$Force.Plate.1_COP.X_mm <- COP(AndrewWong.df, leftorright = "left")$x
AndrewWong.df$Force.Plate.1_COP.Y_mm <- COP(AndrewWong.df, leftorright = "left")$y

AndrewWong.df$Force.Plate.2_COP.X_mm <- COP(AndrewWong.df, leftorright = "right")$x
AndrewWong.df$Force.Plate.2_COP.Y_mm <- COP(AndrewWong.df, leftorright = "right")$y

Daniel.df$Force.Plate.1_COP.X_mm <- COP(Daniel.df, leftorright = "left")$x
Daniel.df$Force.Plate.1_COP.Y_mm <- COP(Daniel.df, leftorright = "left")$y

Daniel.df$Force.Plate.2_COP.X_mm <- COP(Daniel.df, leftorright = "right")$x
Daniel.df$Force.Plate.2_COP.Y_mm <- COP(Daniel.df, leftorright = "right")$y

Frank.df$Force.Plate.1_COP.X_mm <- COP(Frank.df, leftorright = "left")$x
Frank.df$Force.Plate.1_COP.Y_mm <- COP(Frank.df, leftorright = "left")$y

Frank.df$Force.Plate.2_COP.X_mm <- COP(Frank.df, leftorright = "right")$x
Frank.df$Force.Plate.2_COP.Y_mm <- COP(Frank.df, leftorright = "right")$y

Jasmine.df$Force.Plate.1_COP.X_mm <- COP(Jasmine.df, leftorright = "left")$x
Jasmine.df$Force.Plate.1_COP.Y_mm <- COP(Jasmine.df, leftorright = "left")$y

Jasmine.df$Force.Plate.2_COP.X_mm <- COP(Jasmine.df, leftorright = "right")$x
Jasmine.df$Force.Plate.2_COP.Y_mm <- COP(Jasmine.df, leftorright = "right")$y

Jonathan.df$Force.Plate.1_COP.X_mm <- COP(Jonathan.df, leftorright = "left")$x
Jonathan.df$Force.Plate.1_COP.Y_mm <- COP(Jonathan.df, leftorright = "left")$y

Jonathan.df$Force.Plate.2_COP.X_mm <- COP(Jonathan.df, leftorright = "right")$x
Jonathan.df$Force.Plate.2_COP.Y_mm <- COP(Jonathan.df, leftorright = "right")$y

Kelly.df$Force.Plate.1_COP.X_mm <- COP(Kelly.df, leftorright = "left")$x
Kelly.df$Force.Plate.1_COP.Y_mm <- COP(Kelly.df, leftorright = "left")$y

Kelly.df$Force.Plate.2_COP.X_mm <- COP(Kelly.df, leftorright = "right")$x
Kelly.df$Force.Plate.2_COP.Y_mm <- COP(Kelly.df, leftorright = "right")$y

Phoebe.df$Force.Plate.1_COP.X_mm <- COP(Phoebe.df, leftorright = "left")$x
Phoebe.df$Force.Plate.1_COP.Y_mm <- COP(Phoebe.df, leftorright = "left")$y

Phoebe.df$Force.Plate.2_COP.X_mm <- COP(Phoebe.df, leftorright = "right")$x
Phoebe.df$Force.Plate.2_COP.Y_mm <- COP(Phoebe.df, leftorright = "right")$y

Rakesh.df$Force.Plate.1_COP.X_mm <- COP(Rakesh.df, leftorright = "left")$x
Rakesh.df$Force.Plate.1_COP.Y_mm <- COP(Rakesh.df, leftorright = "left")$y

Rakesh.df$Force.Plate.2_COP.X_mm <- COP(Rakesh.df, leftorright = "right")$x
Rakesh.df$Force.Plate.2_COP.Y_mm <- COP(Rakesh.df, leftorright = "right")$y

Thor.df$Force.Plate.1_COP.X_mm <- COP(Thor.df, leftorright = "left")$x
Thor.df$Force.Plate.1_COP.Y_mm <- COP(Thor.df, leftorright = "left")$y

Thor.df$Force.Plate.2_COP.X_mm <- COP(Thor.df, leftorright = "right")$x
Thor.df$Force.Plate.2_COP.Y_mm <- COP(Thor.df, leftorright = "right")$y

Wayne.df$Force.Plate.1_COP.X_mm <- COP(Wayne.df, leftorright = "left")$x
Wayne.df$Force.Plate.1_COP.Y_mm <- COP(Wayne.df, leftorright = "left")$y

Wayne.df$Force.Plate.2_COP.X_mm <- COP(Wayne.df, leftorright = "right")$x
Wayne.df$Force.Plate.2_COP.Y_mm <- COP(Wayne.df, leftorright = "right")$y
#---------------------------------------------------------------------------------------------------

#Calculate FPA--------------------------------------------------------------------------------------
AndrewBollen.df$FPA_Left <- FPA(AndrewBollen.df)$Left
AndrewBollen.df$FPA_Right <- FPA(AndrewBollen.df)$Right

AndrewWong.df$FPA_Left <- FPA(AndrewWong.df)$Left
AndrewWong.df$FPA_Right <- FPA(AndrewWong.df)$Right

Daniel.df$FPA_Left <- FPA(Daniel.df)$Left
Daniel.df$FPA_Right <- FPA(Daniel.df)$Right

Frank.df$FPA_Left <- FPA(Frank.df)$Left
Frank.df$FPA_Right <- FPA(Frank.df)$Right

Jasmine.df$FPA_Left <- FPA(Jasmine.df)$Left
Jasmine.df$FPA_Right <- FPA(Jasmine.df)$Right

Jonathan.df$FPA_Left <- FPA(Jonathan.df)$Left
Jonathan.df$FPA_Right <- FPA(Jonathan.df)$Right

Kelly.df$FPA_Left <- FPA(Kelly.df)$Left
Kelly.df$FPA_Right <- FPA(Kelly.df)$Right

Phoebe.df$FPA_Left <- FPA(Phoebe.df)$Left
Phoebe.df$FPA_Right <- FPA(Phoebe.df)$Right

Rakesh.df$FPA_Left <- FPA(Rakesh.df)$Left
Rakesh.df$FPA_Right <- FPA(Rakesh.df)$Right

Thor.df$FPA_Left <- FPA(Thor.df)$Left
Thor.df$FPA_Right <- FPA(Thor.df)$Right

Wayne.df$FPA_Left <- FPA(Wayne.df)$Left
Wayne.df$FPA_Right <- FPA(Wayne.df)$Right
#---------------------------------------------------------------------------------------------------

#Calculate Calcaneal Marker during Initial Step-----------------------------------------------------
AndrewBollen.df$StartLeftStep <- Step(AndrewBollen.df)$StartLeft
AndrewBollen.df$StartRightStep <- Step(AndrewBollen.df)$StartRight
AndrewBollen.df$EndLeftStep <- Step(AndrewBollen.df)$EndLeft
AndrewBollen.df$EndRightStep <- Step(AndrewBollen.df)$EndRight

AndrewWong.df$StartLeftStep <- Step(AndrewWong.df)$StartLeft
AndrewWong.df$StartRightStep <- Step(AndrewWong.df)$StartRight
AndrewWong.df$EndLeftStep <- Step(AndrewWong.df)$EndLeft
AndrewWong.df$EndRightStep <- Step(AndrewWong.df)$EndRight

Daniel.df$StartLeftStep <- Step(Daniel.df)$StartLeft
Daniel.df$StartRightStep <- Step(Daniel.df)$StartRight
Daniel.df$EndLeftStep <- Step(Daniel.df)$EndLeft
Daniel.df$EndRightStep <- Step(Daniel.df)$EndRight

Frank.df$StartLeftStep <- Step(Frank.df)$StartLeft
Frank.df$StartRightStep <- Step(Frank.df)$StartRight
Frank.df$EndLeftStep <- Step(Frank.df)$EndLeft
Frank.df$EndRightStep <- Step(Frank.df)$EndRight

Jasmine.df$StartLeftStep <- Step(Jasmine.df)$StartLeft
Jasmine.df$StartRightStep <- Step(Jasmine.df)$StartRight
Jasmine.df$EndLeftStep <- Step(Jasmine.df)$EndLeft
Jasmine.df$EndRightStep <- Step(Jasmine.df)$EndRight

Jonathan.df$StartLeftStep <- Step(Jonathan.df)$StartLeft
Jonathan.df$StartRightStep <- Step(Jonathan.df)$StartRight
Jonathan.df$EndLeftStep <- Step(Jonathan.df)$EndLeft
Jonathan.df$EndRightStep <- Step(Jonathan.df)$EndRight

Kelly.df$StartLeftStep <- Step(Kelly.df)$StartLeft
Kelly.df$StartRightStep <- Step(Kelly.df)$StartRight
Kelly.df$EndLeftStep <- Step(Kelly.df)$EndLeft
Kelly.df$EndRightStep <- Step(Kelly.df)$EndRight

Phoebe.df$StartLeftStep <- Step(Phoebe.df)$StartLeft
Phoebe.df$StartRightStep <- Step(Phoebe.df)$StartRight
Phoebe.df$EndLeftStep <- Step(Phoebe.df)$EndLeft
Phoebe.df$EndRightStep <- Step(Phoebe.df)$EndRight

Rakesh.df$StartLeftStep <- Step(Rakesh.df)$StartLeft
Rakesh.df$StartRightStep <- Step(Rakesh.df)$StartRight
Rakesh.df$EndLeftStep <- Step(Rakesh.df)$EndLeft
Rakesh.df$EndRightStep <- Step(Rakesh.df)$EndRight

Thor.df$StartLeftStep <- Step(Thor.df)$StartLeft
Thor.df$StartRightStep <- Step(Thor.df)$StartRight
Thor.df$EndLeftStep <- Step(Thor.df)$EndLeft
Thor.df$EndRightStep <- Step(Thor.df)$EndRight

Wayne.df$StartLeftStep <- Step(Wayne.df)$StartLeft
Wayne.df$StartRightStep <- Step(Wayne.df)$StartRight
Wayne.df$EndLeftStep <- Step(Wayne.df)$EndLeft
Wayne.df$EndRightStep <- Step(Wayne.df)$EndRight
#---------------------------------------------------------------------------------------------------

#Calculate Step Width SW----------------------------------------------------------------------------
AndrewBollen.df$SW <- StepWidth(AndrewBollen.df)

AndrewWong.df$SW <- StepWidth(AndrewWong.df)

Daniel.df$SW <- StepWidth(Daniel.df)

Frank.df$SW <- StepWidth(Frank.df)

Jasmine.df$SW <- StepWidth(Jasmine.df)

Jonathan.df$SW <- StepWidth(Jonathan.df)

Kelly.df$SW <- StepWidth(Kelly.df)

Phoebe.df$SW <- StepWidth(Phoebe.df)

Rakesh.df$SW <- StepWidth(Rakesh.df)

Thor.df$SW <- StepWidth(Thor.df)

Wayne.df$SW <- StepWidth(Wayne.df)
#---------------------------------------------------------------------------------------------------

#Calculate position vector from knee midpoint to COP-------------------------------------------------------------------------------------
AndrewBollen.df$LposX <- PosVec(AndrewBollen.df)$Lx
AndrewBollen.df$LposY <- PosVec(AndrewBollen.df)$Ly
AndrewBollen.df$LposZ <- PosVec(AndrewBollen.df)$Lz
AndrewBollen.df$RposX <- PosVec(AndrewBollen.df)$Rx
AndrewBollen.df$RposY <- PosVec(AndrewBollen.df)$Ry
AndrewBollen.df$RposZ <- PosVec(AndrewBollen.df)$Rz

AndrewWong.df$LposX <- PosVec(AndrewWong.df)$Lx
AndrewWong.df$LposY <- PosVec(AndrewWong.df)$Ly
AndrewWong.df$LposZ <- PosVec(AndrewWong.df)$Lz
AndrewWong.df$RposX <- PosVec(AndrewWong.df)$Rx
AndrewWong.df$RposY <- PosVec(AndrewWong.df)$Ry
AndrewWong.df$RposZ <- PosVec(AndrewWong.df)$Rz

Daniel.df$LposX <- PosVec(Daniel.df)$Lx
Daniel.df$LposY <- PosVec(Daniel.df)$Ly
Daniel.df$LposZ <- PosVec(Daniel.df)$Lz
Daniel.df$RposX <- PosVec(Daniel.df)$Rx
Daniel.df$RposY <- PosVec(Daniel.df)$Ry
Daniel.df$RposZ <- PosVec(Daniel.df)$Rz

Frank.df$LposX <- PosVec(Frank.df)$Lx
Frank.df$LposY <- PosVec(Frank.df)$Ly
Frank.df$LposZ <- PosVec(Frank.df)$Lz
Frank.df$RposX <- PosVec(Frank.df)$Rx
Frank.df$RposY <- PosVec(Frank.df)$Ry
Frank.df$RposZ <- PosVec(Frank.df)$Rz

Jasmine.df$LposX <- PosVec(Jasmine.df)$Lx
Jasmine.df$LposY <- PosVec(Jasmine.df)$Ly
Jasmine.df$LposZ <- PosVec(Jasmine.df)$Lz
Jasmine.df$RposX <- PosVec(Jasmine.df)$Rx
Jasmine.df$RposY <- PosVec(Jasmine.df)$Ry
Jasmine.df$RposZ <- PosVec(Jasmine.df)$Rz

Jonathan.df$LposX <- PosVec(Jonathan.df)$Lx
Jonathan.df$LposY <- PosVec(Jonathan.df)$Ly
Jonathan.df$LposZ <- PosVec(Jonathan.df)$Lz
Jonathan.df$RposX <- PosVec(Jonathan.df)$Rx
Jonathan.df$RposY <- PosVec(Jonathan.df)$Ry
Jonathan.df$RposZ <- PosVec(Jonathan.df)$Rz

Kelly.df$LposX <- PosVec(Kelly.df)$Lx
Kelly.df$LposY <- PosVec(Kelly.df)$Ly
Kelly.df$LposZ <- PosVec(Kelly.df)$Lz
Kelly.df$RposX <- PosVec(Kelly.df)$Rx
Kelly.df$RposY <- PosVec(Kelly.df)$Ry
Kelly.df$RposZ <- PosVec(Kelly.df)$Rz

Phoebe.df$LposX <- PosVec(Phoebe.df)$Lx
Phoebe.df$LposY <- PosVec(Phoebe.df)$Ly
Phoebe.df$LposZ <- PosVec(Phoebe.df)$Lz
Phoebe.df$RposX <- PosVec(Phoebe.df)$Rx
Phoebe.df$RposY <- PosVec(Phoebe.df)$Ry
Phoebe.df$RposZ <- PosVec(Phoebe.df)$Rz

Rakesh.df$LposX <- PosVec(Rakesh.df)$Lx
Rakesh.df$LposY <- PosVec(Rakesh.df)$Ly
Rakesh.df$LposZ <- PosVec(Rakesh.df)$Lz
Rakesh.df$RposX <- PosVec(Rakesh.df)$Rx
Rakesh.df$RposY <- PosVec(Rakesh.df)$Ry
Rakesh.df$RposZ <- PosVec(Rakesh.df)$Rz

Thor.df$LposX <- PosVec(Thor.df)$Lx
Thor.df$LposY <- PosVec(Thor.df)$Ly
Thor.df$LposZ <- PosVec(Thor.df)$Lz
Thor.df$RposX <- PosVec(Thor.df)$Rx
Thor.df$RposY <- PosVec(Thor.df)$Ry
Thor.df$RposZ <- PosVec(Thor.df)$Rz

Wayne.df$LposX <- PosVec(Wayne.df)$Lx
Wayne.df$LposY <- PosVec(Wayne.df)$Ly
Wayne.df$LposZ <- PosVec(Wayne.df)$Lz
Wayne.df$RposX <- PosVec(Wayne.df)$Rx
Wayne.df$RposY <- PosVec(Wayne.df)$Ry
Wayne.df$RposZ <- PosVec(Wayne.df)$Rz
#----------------------------------------------------------------------------------------------------------------------------------------

#Calculate KAM--------------------------------------------------------
AndrewBollen.df$LKAM <- CalculateKam(AndrewBollen.df)$Left
AndrewBollen.df$RKAM <- CalculateKam(AndrewBollen.df)$Right

AndrewWong.df$LKAM <- CalculateKam(AndrewWong.df)$Left
AndrewWong.df$RKAM <- CalculateKam(AndrewWong.df)$Right

Daniel.df$LKAM <- CalculateKam(Daniel.df)$Left
Daniel.df$RKAM <- CalculateKam(Daniel.df)$Right

Frank.df$LKAM <- CalculateKam(Frank.df)$Left
Frank.df$RKAM <- CalculateKam(Frank.df)$Right

Jasmine.df$LKAM <- CalculateKam(Jasmine.df)$Left
Jasmine.df$RKAM <- CalculateKam(Jasmine.df)$Right

Jonathan.df$LKAM <- CalculateKam(Jonathan.df)$Left
Jonathan.df$RKAM <- CalculateKam(Jonathan.df)$Right

Kelly.df$LKAM <- CalculateKam(Kelly.df)$Left
Kelly.df$RKAM <- CalculateKam(Kelly.df)$Right

Phoebe.df$LKAM <- CalculateKam(Phoebe.df)$Left
Phoebe.df$RKAM <- CalculateKam(Phoebe.df)$Right

Rakesh.df$LKAM <- CalculateKam(Rakesh.df)$Left
Rakesh.df$RKAM <- CalculateKam(Rakesh.df)$Right

Thor.df$LKAM <- CalculateKam(Thor.df)$Left
Thor.df$RKAM <- CalculateKam(Thor.df)$Right

Wayne.df$LKAM <- CalculateKam(Wayne.df)$Left
Wayne.df$RKAM <- CalculateKam(Wayne.df)$Right
#----------------------------------------------------------------------------------------------------------------------------------------

#Replace NAs with 0-----------------------------------------
AndrewBollen.df[is.na(AndrewBollen.df)] <- 0
AndrewWong.df[is.na(AndrewWong.df)] <- 0
Daniel.df[is.na(Daniel.df)] <- 0
Frank.df[is.na(Frank.df)] <- 0
Jasmine.df[is.na(Jasmine.df)] <- 0
Jonathan.df[is.na(Jonathan.df)] <- 0
Kelly.df[is.na(Kelly.df)] <- 0
Phoebe.df[is.na(Phoebe.df)] <- 0
Rakesh.df[is.na(Rakesh.df)] <- 0
Thor.df[is.na(Thor.df)] <- 0
Wayne.df[is.na(Wayne.df)] <- 0
#-----------------------------------------------------------

#Normalize KAM in Nm with BW-------------------------------------------------------------------------------------------------------------------
AndrewBollen.df$LKAM_norm <- (AndrewBollen.df$LKAM/1000)/(72*1.72)
AndrewBollen.df$RKAM_norm <- (AndrewBollen.df$RKAM/1000)/(72*1.72)

AndrewWong.df$LKAM_norm <- (AndrewWong.df$LKAM/1000)/(66*1.76)
AndrewWong.df$RKAM_norm <- (AndrewWong.df$RKAM/1000)/(66*1.76)

Daniel.df$LKAM_norm <- (Daniel.df$LKAM/1000)/(73*1.74)
Daniel.df$RKAM_norm <- (Daniel.df$RKAM/1000)/(73*1.74)

Frank.df$LKAM_norm <- (Frank.df$LKAM/1000)/(71.5*1.74)
Frank.df$RKAM_norm <- (Frank.df$RKAM/1000)/(71.5*1.74)

Jasmine.df$LKAM_norm <- (Jasmine.df$LKAM/1000)/(58*1.67)
Jasmine.df$RKAM_norm <- (Jasmine.df$RKAM/1000)/(58*1.67)

Jonathan.df$LKAM_norm <- (Jonathan.df$LKAM/1000)/(75*1.85)
Jonathan.df$RKAM_norm <- (Jonathan.df$RKAM/1000)/(75*1.85)

Kelly.df$LKAM_norm <- (Kelly.df$LKAM/1000)/(77*1.69)
Kelly.df$RKAM_norm <- (Kelly.df$RKAM/1000)/(77*1.69)

Phoebe.df$LKAM_norm <- (Phoebe.df$LKAM/1000)/(51*1.6)
Phoebe.df$RKAM_norm <- (Phoebe.df$RKAM/1000)/(51*1.6)

Rakesh.df$LKAM_norm <- (Rakesh.df$LKAM/1000)/(91*1.75)
Rakesh.df$RKAM_norm <- (Rakesh.df$RKAM/1000)/(91*1.75)

Thor.df$LKAM_norm <- (Thor.df$LKAM/1000)/(74*1.76) 
Thor.df$RKAM_norm <- (Thor.df$RKAM/1000)/(74*1.76)

Wayne.df$LKAM_norm <- (Wayne.df$LKAM/1000)/(95*1.70)
Wayne.df$RKAM_norm <- (Wayne.df$RKAM/1000)/(95*1.70)
#----------------------------------------------------------------------------------------------------------------------------------------

#Remove data from dataframe which is no longer needed-------------------------------------
dataframeColumns <- c('Force.Plate.1_Force.Z_N', 'Force.Plate.2_Force.Z_N', 'Force.Plate.1_COP.X_mm', 'Force.Plate.1_COP.Y_mm', 'Force.Plate.2_COP.X_mm', 'Force.Plate.2_COP.Y_mm', 'FPA_Left', 'FPA_Right', 'SW', 'StartLeftStep', 'StartRightStep', 'EndLeftStep', 'EndRightStep', 'LKAM_norm', 'RKAM_norm', 'LKAM', 'RKAM')

AndrewBollen.df <- AndrewBollen.df[,dataframeColumns]
AndrewWong.df <- AndrewWong.df[,dataframeColumns]
Daniel.df <- Daniel.df[,dataframeColumns]
Frank.df <- Frank.df[,dataframeColumns]
Jasmine.df <- Jasmine.df[,dataframeColumns]
Jonathan.df <- Jonathan.df[,dataframeColumns]
Kelly.df <- Kelly.df[,dataframeColumns]
Phoebe.df <- Phoebe.df[,dataframeColumns]
Rakesh.df <- Rakesh.df[,dataframeColumns]
Thor.df <- Thor.df[,dataframeColumns]
Wayne.df <- Wayne.df[,dataframeColumns]
#-----------------------------------------------------------------------------------------

#Check KAM with Start and End step and check where they need clipping---------------------
#Plot line graph of KAM and Start and End
Plot <- plotSteps(AndrewBollen.df, leftorright = "left")
Plot
Plot <- plotSteps(AndrewBollen.df, leftorright = "right")
Plot

Plot <- plotSteps(AndrewWong.df, leftorright = "left")
Plot
Plot <- plotSteps(AndrewWong.df, leftorright = "right")
Plot

Plot <- plotSteps(Daniel.df, leftorright = "left")
Plot
Plot <- plotSteps(Daniel.df, leftorright = "right")
Plot

Plot <- plotSteps(Frank.df, leftorright = "left")
Plot
Plot <- plotSteps(Frank.df, leftorright = "right")
Plot

Plot <- plotSteps(Jasmine.df, leftorright = "left")
Plot
Plot <- plotSteps(Jasmine.df, leftorright = "right")
Plot

Plot <- plotSteps(Jonathan.df, leftorright = "left")
Plot
Plot <- plotSteps(Jonathan.df, leftorright = "right")
Plot

Plot <- plotSteps(Kelly.df, leftorright = "left")
Plot
Plot <- plotSteps(Kelly.df, leftorright = "right")
Plot

Plot <- plotSteps(Phoebe.df, leftorright = "left")
Plot
Plot <- plotSteps(Phoebe.df, leftorright = "right")
Plot

Plot <- plotSteps(Rakesh.df, leftorright = "left")
Plot
Plot <- plotSteps(Rakesh.df, leftorright = "right")
Plot

Plot <- plotSteps(Thor.df, leftorright = "left")
Plot
Plot <- plotSteps(Thor.df, leftorright = "right")
Plot

Plot <- plotSteps(Wayne.df, leftorright = "left")
Plot
Plot <- plotSteps(Wayne.df, leftorright = "right")
Plot
#-----------------------------------------------------------------------------------------

#Clip dataframes to just relevant sections------------------------------------------------
AndrewBollen.df <- clipAt(AndrewBollen.df, leftorright = "left", 3, 2)
AndrewBollen.df <- clipAt(AndrewBollen.df, leftorright = "right", 3, 3)

AndrewWong.df <- clipAt(AndrewWong.df, leftorright = "left", 2, 1)
AndrewWong.df <- clipAt(AndrewWong.df, leftorright = "right", 4, 1)

Daniel.df <- clipAt(Daniel.df, leftorright = "left", 2, 5)
Daniel.df <- clipAt(Daniel.df, leftorright = "right", 2, 5)

Frank.df <- clipAt(Frank.df, leftorright = "left", 2, 1)
Frank.df <- clipAt(Frank.df, leftorright = "right", 2, 2)

Jasmine.df <- clipAt(Jasmine.df, leftorright = "left", 3, 3)
Jasmine.df <- clipAt(Jasmine.df, leftorright = "right", 2, 1)

Jonathan.df <- clipAt(Jonathan.df, leftorright = "left", 2, 3)
Jonathan.df <- clipAt(Jonathan.df, leftorright = "right", 3, 3)

Kelly.df <- clipAt(Kelly.df, leftorright = "left", 2, 15)
Kelly.df <- clipAt(Kelly.df, leftorright = "right", 2, 14)

Phoebe.df <- clipAt(Phoebe.df, leftorright = "left", 3, 2)
Phoebe.df <- clipAt(Phoebe.df, leftorright = "right", 4, 2)

Rakesh.df <- clipAt(Rakesh.df, leftorright = "left", 2, 1)
Rakesh.df <- clipAt(Rakesh.df, leftorright = "right", 2, 1)

Thor.df <- clipAt(Thor.df, leftorright = "left", 2, 1)
Thor.df <- clipAt(Thor.df, leftorright = "right", 2, 1)

Wayne.df <- clipAt(Wayne.df, leftorright = "left", 2, 1)
Wayne.df <- clipAt(Wayne.df, leftorright = "right", 2, 1)
#-----------------------------------------------------------------------------------------

#Identify first and second peaks with their indices---------------------------------------
AndrewBollenLeftPeak <- PeakSelect(AndrewBollen.df, leftorright = "left") 
AndrewBollenRightPeak <- PeakSelect(AndrewBollen.df, leftorright = "right")

AndrewWongLeftPeak <- PeakSelect(AndrewWong.df, leftorright = "left")
AndrewWongRightPeak <- PeakSelect(AndrewWong.df, leftorright = "right")

DanielLeftPeak <- PeakSelect(Daniel.df, leftorright = "left")
DanielRightPeak <- PeakSelect(Daniel.df, leftorright = "right")

FrankLeftPeak <- PeakSelect(Frank.df, leftorright = "left")
FrankRightPeak <- PeakSelect(Frank.df, leftorright = "right")

JasmineLeftPeak <- PeakSelect(Jasmine.df, leftorright = "left")
JasmineRightPeak <- PeakSelect(Jasmine.df, leftorright = "right")

JonathanLeftPeak <- PeakSelect(Jonathan.df, leftorright = "left")
JonathanRightPeak <- PeakSelect(Jonathan.df, leftorright = "right")

KellyLeftPeak <- PeakSelect(Kelly.df, leftorright = "left")
KellyRightPeak <- PeakSelect(Kelly.df, leftorright = "right")

PhoebeLeftPeak <- PeakSelect(Phoebe.df, leftorright = "left")
PhoebeRightPeak <- PeakSelect(Phoebe.df, leftorright = "right")

RakeshLeftPeak <- PeakSelect(Rakesh.df, leftorright = "left")
RakeshRightPeak <- PeakSelect(Rakesh.df, leftorright = "right")

ThorLeftPeak <- PeakSelect(Thor.df, leftorright = "left")
ThorRightPeak <- PeakSelect(Thor.df, leftorright = "right")

WayneLeftPeak <- PeakSelect(Wayne.df, leftorright = "left")
WayneRightPeak <- PeakSelect(Wayne.df, leftorright = "right")
#-----------------------------------------------------------------------------------------

dir <- "C:\\Users\\bioengsu\\Desktop\\KAM\\Subject KAM CSV"
setwd(dir)
#Create dataframe with Subject's FPA, SW, first and second peak KAM at every step---------
AndrewBollenLeftStep <- peakDataframe(AndrewBollen.df, AndrewBollenLeftPeak,leftorright = "left")
AndrewBollenRightStep <- peakDataframe(AndrewBollen.df, AndrewBollenRightPeak,leftorright = "right")
write.csv(AndrewBollenLeftStep, file = "AndrewBollenLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(AndrewBollenRightStep, file = "AndrewBollenRightStep.csv", row.names = FALSE, sep = ",")

AndrewWongLeftStep <- peakDataframe(AndrewWong.df, AndrewWongLeftPeak,leftorright = "left")
AndrewWongRightStep <- peakDataframe(AndrewWong.df, AndrewWongRightPeak,leftorright = "right")
write.csv(AndrewWongLeftStep, file = "AndrewWongLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(AndrewWongRightStep, file = "AndrewWongRightStep.csv", row.names = FALSE, sep = ",")

DanielLeftStep <- peakDataframe(Daniel.df, DanielLeftPeak,leftorright = "left")
DanielRightStep <- peakDataframe(Daniel.df, DanielRightPeak,leftorright = "right")
write.csv(DanielLeftStep, file = "DanielLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(DanielRightStep, file = "DanielRightStep.csv", row.names = FALSE, sep = ",")

FrankLeftStep <- peakDataframe(Frank.df, FrankLeftPeak,leftorright = "left")
FrankRightStep <- peakDataframe(Frank.df, FrankRightPeak,leftorright = "right")
write.csv(FrankLeftStep, file = "FrankLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(FrankRightStep, file = "FrankRightStep.csv", row.names = FALSE, sep = ",")

JasmineLeftStep <- peakDataframe(Jasmine.df, JasmineLeftPeak,leftorright = "left")
JasmineRightStep <- peakDataframe(Jasmine.df, JasmineRightPeak,leftorright = "right")
write.csv(JasmineLeftStep, file = "JasmineLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(JasmineRightStep, file = "JasmineRightStep.csv", row.names = FALSE, sep = ",")

JonathanLeftStep <- peakDataframe(Jonathan.df, JonathanLeftPeak,leftorright = "left")
JonathanRightStep <- peakDataframe(Jonathan.df, JonathanRightPeak,leftorright = "right")
write.csv(JonathanLeftStep, file = "JonathanLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(JonathanRightStep, file = "JonathanRightStep.csv", row.names = FALSE, sep = ",")

KellyLeftStep <- peakDataframe(Kelly.df, KellyLeftPeak,leftorright = "left")
KellyRightStep <- peakDataframe(Kelly.df, KellyRightPeak,leftorright = "right")
write.csv(KellyLeftStep, file = "KellyLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(KellyRightStep, file = "KellyRightStep.csv", row.names = FALSE, sep = ",")

PhoebeLeftStep <- peakDataframe(Phoebe.df, PhoebeLeftPeak,leftorright = "left")
PhoebeRightStep <- peakDataframe(Phoebe.df, PhoebeRightPeak,leftorright = "right")
write.csv(PhoebeLeftStep, file = "PhoebeLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(PhoebeRightStep, file = "PhoebeRightStep.csv", row.names = FALSE, sep = ",")

RakeshLeftStep <- peakDataframe(Rakesh.df, RakeshLeftPeak,leftorright = "left")
RakeshRightStep <- peakDataframe(Rakesh.df, RakeshRightPeak,leftorright = "right")
write.csv(RakeshLeftStep, file = "RakeshLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(RakeshRightStep, file = "RakeshRightStep.csv", row.names = FALSE, sep = ",")

ThorLeftStep <- peakDataframe(Thor.df, ThorLeftPeak,leftorright = "left")
ThorRightStep <- peakDataframe(Thor.df, ThorRightPeak,leftorright = "right")
write.csv(ThorLeftStep, file = "ThorLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(ThorRightStep, file = "ThorRightStep.csv", row.names = FALSE, sep = ",")

WayneLeftStep <- peakDataframe(Wayne.df, WayneLeftPeak,leftorright = "left")
WayneRightStep <- peakDataframe(Wayne.df, WayneRightPeak,leftorright = "right")
write.csv(WayneLeftStep, file = "WayneLeftStep.csv", row.names = FALSE, sep = ",")
write.csv(WayneRightStep, file = "WayneRightStep.csv", row.names = FALSE, sep = ",")
#-----------------------------------------------------------------------------------------