#updated KAM script for Daniel and Thor

rm(list = ls())

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

#functions-----------------------------------------------
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

makeZero <- function(dataframe){
  attach(dataframe)
  
  for (i in 1:nrow(dataframe)){
    if (Force.Plate.1_Force.Z_N[i] < 80){
      dataframe[i,25:36] <- 0
    }
    if (Force.Plate.2_Force.Z_N[i] < 80){
      dataframe[i,36:ncol(dataframe)] <- 0
    }
  }
  
  detach(dataframe)
  
  return(dataframe)
}

#Calculates COP from force plate data
#(238.5, 812.5, 0) for force plate 1, (727.5, 812.5, 0) for force plate 2??
COP <- function(dataframe, leftorright = "left"){
  attach(dataframe)
  h <- 15 #height of padding
  
  if (leftorright == "left"){ 
    COPx <- ((-h*Force.Plate.1_Force.X_N - Force.Plate.1_Moment.Y_Nmm)/Force.Plate.1_Force.Z_N) + 238.5
    COPy <- ((h*Force.Plate.1_Force.Y_N + Force.Plate.1_Moment.X_Nmm)/Force.Plate.1_Force.Z_N) + 812.5
    
  } else if (leftorright == "right"){
    COPx <- ((-h*Force.Plate.2_Force.X_N - Force.Plate.2_Moment.Y_Nmm)/Force.Plate.2_Force.Z_N) + 727.5
    COPy <- ((h*Force.Plate.2_Force.Y_N + Force.Plate.2_Moment.X_Nmm)/Force.Plate.2_Force.Z_N) + 812.5
  }
  
  COPz <- 0

  detach(dataframe)
  COP <- list("x" = COPx, "y" = COPy, "z" = COPz)
  
  return(COP)
}

filterDataframe <- function(dataframe){
  
  #Replace NAs with 0
  dataframe[is.na(dataframe)] <- 0
  
  for (i in 1:ncol(dataframe)){
    dataframe[,i] <- with(data = dataframe, expr = filterSignal(dataframe[,i], cutoff = 10, type = 'low', sampleFreq = 100))
  }
  
  return(dataframe)
}

#Filter signal
filterSignal <- function(signal, cutoff, type, sampleFreq){
  
  bf <- butter(4, (2*cutoff)/(sampleFreq), type = type)
  y <- filtfilt(bf, signal)
  
  return(y)
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

#Calculate and update dataframe with FPA
FPA <- function(dataframe){
  LeftDot <- c()
  RightDot <- c()
  
  #Calculate Calcaneal to Toe vector
  LeftFootVec_X <- positionVec(origin = dataframe$LC_X, endPoint = dataframe$LT_X)
  RightFootVec_X <- positionVec(origin = dataframe$RC_X, endPoint = dataframe$RT_X)
  
  LeftFootVec_Y <- positionVec(origin = dataframe$LC_Y, endPoint = dataframe$LT_Y)
  RightFootVec_Y <- positionVec(origin = dataframe$RC_Y, endPoint = dataframe$RT_Y)
  
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
      LeftHeel[i] <- dataframe[i,'LC_X']
    } else {
      if (i != 1){
        LeftHeel[i] <- LeftHeel[i-1]
      }
    }
    
    if (dataframe[i,'StartRightStep']){
      RightHeel[i] <- dataframe[i,'RC_X']
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

CalculateKam <- function(dataframe, type = "actual"){
  
  headers <- colnames(dataframe)
  
  # for (i in 1:length(headers)){
  #   if ((headers[i] == "LTibiaToGlobal_Z_mm")||(headers[i] == "LTibiaToGlobal_Z_degrees")||(headers[i] == "LTibiaToGlobal_Z")){
  #     headers[i] <- "LTibiaToGlobal_Z_deg"
  #   }
  #   
  #   if ((headers[i] == "RTibiaToGlobal_Z_mm")||(headers[i] == "RTibiaToGlobal_Z_degrees")||(headers[i] == "RTibiaToGlobal_Z")){
  #     headers[i] <- "RTibiaToGlobal_Z_deg"
  #   }
  # }
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
      if (type == "actual"){
        KAMx_Left <- append(KAMx_Left, KAM(c(LposX[i], LposY[i], LposZ[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[1])
        KAMy_Left <- append(KAMy_Left, KAM(c(LposX[i], LposY[i], LposZ[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[2])
        KAMz_Left <- append(KAMz_Left, KAM(c(LposX[i], LposY[i], LposZ[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[3])
        
        #Rotate KAM to Tibia Frame
        LTibKAMx <- append(LTibKAMx, rotate(LTibiaToGlobal_Z[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[1])
        LTibKAMy <- append(LTibKAMy, rotate(LTibiaToGlobal_Z[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[2])
        LTibKAMz <- append(LTibKAMz, rotate(LTibiaToGlobal_Z[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[3])
      } else if (type == "vicon"){
        KAMx_Left <- append(KAMx_Left, KAM(c(LposX_V[i], LposY_V[i], LposZ_V[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[1])
        KAMy_Left <- append(KAMy_Left, KAM(c(LposX_V[i], LposY_V[i], LposZ_V[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[2])
        KAMz_Left <- append(KAMz_Left, KAM(c(LposX_V[i], LposY_V[i], LposZ_V[i]), c(Force.Plate.1_Force.X_N[i], Force.Plate.1_Force.Y_N[i], Force.Plate.1_Force.Z_N[i]))[3])
        
        #Rotate KAM to Tibia Frame
        LTibKAMx <- append(LTibKAMx, rotate(LTibiaToGlobal_Z[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[1])
        LTibKAMy <- append(LTibKAMy, rotate(LTibiaToGlobal_Z[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[2])
        LTibKAMz <- append(LTibKAMz, rotate(LTibiaToGlobal_Z[i], KAMx_Left[i], KAMy_Left[i], KAMz_Left[i])[3])
      }
    } else {
      LTibKAMx <- append(LTibKAMx, 0)
      LTibKAMy <- append(LTibKAMy, 0)
      LTibKAMz <- append(LTibKAMz, 0)
      
      KAMx_Left <- append(KAMx_Left, 0)
      KAMy_Left <- append(KAMy_Left, 0)
      KAMz_Left <- append(KAMz_Left, 0)
    }
    
    if (Force.Plate.2_Force.Z_N[i] > 80){
      if (type == "actual"){
        KAMx_Right <- append(KAMx_Right, KAM(c(RposX[i], RposY[i], RposZ[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[1])
        KAMy_Right <- append(KAMy_Right, KAM(c(RposX[i], RposY[i], RposZ[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[2])
        KAMz_Right <- append(KAMz_Right, KAM(c(RposX[i], RposY[i], RposZ[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[3])
        
        #Rotate KAM to Tibia Frame
        RTibKAMx <- append(RTibKAMx, rotate(RTibiaToGlobal_Z[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[1])
        RTibKAMy <- append(RTibKAMy, rotate(RTibiaToGlobal_Z[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[2])
        RTibKAMz <- append(RTibKAMz, rotate(RTibiaToGlobal_Z[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[3])
      } else if (type == "vicon"){
        KAMx_Right <- append(KAMx_Right, KAM(c(RposX_V[i], RposY_V[i], RposZ_V[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[1])
        KAMy_Right <- append(KAMy_Right, KAM(c(RposX_V[i], RposY_V[i], RposZ_V[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[2])
        KAMz_Right <- append(KAMz_Right, KAM(c(RposX_V[i], RposY_V[i], RposZ_V[i]), c(Force.Plate.2_Force.X_N[i], Force.Plate.2_Force.Y_N[i], Force.Plate.2_Force.Z_N[i]))[3])
        
        #Rotate KAM to Tibia Frame
        RTibKAMx <- append(RTibKAMx, rotate(RTibiaToGlobal_Z[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[1])
        RTibKAMy <- append(RTibKAMy, rotate(RTibiaToGlobal_Z[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[2])
        RTibKAMz <- append(RTibKAMz, rotate(RTibiaToGlobal_Z[i], KAMx_Right[i], KAMy_Right[i], KAMz_Right[i])[3])
      }
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
  
  LposX <- mapply(positionVec, dataframe$LKJC_X, dataframe$Force.Plate.1_COP.X_mm)
  LposY <- mapply(positionVec, dataframe$LKJC_Y, dataframe$Force.Plate.1_COP.Y_mm)
  LposZ <- mapply(positionVec, dataframe$LKJC_Z, dataframe$Force.Plate.1_COP.Z_mm)
  
  RposX <- mapply(positionVec, dataframe$RKJC_X, dataframe$Force.Plate.2_COP.X_mm)
  RposY <- mapply(positionVec, dataframe$RKJC_Y, dataframe$Force.Plate.2_COP.Y_mm)
  RposZ <- mapply(positionVec, dataframe$RKJC_Z, dataframe$Force.Plate.2_COP.Z_mm)
  
  LposX_Vicon <- mapply(positionVec, dataframe$LKJC_X, dataframe$GRW1_mm_X)
  LposY_Vicon <- mapply(positionVec, dataframe$LKJC_Y, dataframe$GRW1_mm_Y)
  LposZ_Vicon <- mapply(positionVec, dataframe$LKJC_Z, dataframe$Force.Plate.1_COP.Z_mm)
  
  RposX_Vicon <- mapply(positionVec, dataframe$RKJC_X, dataframe$GRW2_mm_X)
  RposY_Vicon <- mapply(positionVec, dataframe$RKJC_Y, dataframe$GRW2_mm_Y)
  RposZ_Vicon <- mapply(positionVec, dataframe$RKJC_Z, dataframe$Force.Plate.2_COP.Z_mm)
  
  PosList <- list("Lx" = LposX, "Ly" = LposY, "Lz" = LposZ, "Rx" = RposX, "Ry" = RposY, "Rz" = RposZ, 
                  "Lx_V" = LposX_Vicon, "Ly_V" = LposY_Vicon, "Lz_V" = LposZ_Vicon, 
                  "Rx_V" = RposX_Vicon, "Ry_V" = RposY_Vicon, "Rz_V" = RposZ_Vicon)
  
  return(PosList)
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

PeakSelect <- function(dataframe, leftorright = "left"){
  attach(dataframe)
  peaks <- c()
  startIndex <- 0
  endIndex <- 0
  firstPeak <- c()
  secondPeak <- c()
  par(mfrow=c(1,1))
  
  if (leftorright == "left"){
    for (i in 1:nrow(dataframe)){
      if (StartLeftStep[i]){
        startIndex <- i
      }
      
      if (EndLeftStep[i] && (startIndex != 0)){
        endIndex <- i
        print(plot(x = startIndex:endIndex, y = LKAM_norm[startIndex:endIndex], type = 'l'))
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
      
      if (EndRightStep[i] && (startIndex != 0)){
        endIndex <- i
        print(plot(x = startIndex:endIndex, y = RKAM_norm[startIndex:endIndex], type = 'l'))
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
  
  detach(dataframe)
  return(peakers)
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
#--------------------------------------------------------

#read cleaned csv----------------------------------------
dir <- c("C:\\Users\\bioengsu\\Desktop\\KAM\\Cleaned CSV\\Trial 2\\Daniel", "C:\\Users\\bioengsu\\Desktop\\KAM\\Cleaned CSV\\Trial 2\\Thor")

setwd(dir[1])

DanielFiles <- c(choose.files())

DanielCSV <- lapply(DanielFiles, read.csv, header = TRUE, stringsAsFactors = FALSE)

setwd(dir[2])

ThorFiles <- c(choose.files())

ThorCSV <- lapply(ThorFiles, read.csv, header = TRUE, stringsAsFactors = FALSE)
#--------------------------------------------------------

#Keep only Variables of interest------------------------------------------------------------------------------
#RC, RT, LC, LT, RTIBTOGLOBAL_Z, LTIBTOGLOBAL_Z, RKJC, LKJC, FORCE, MOMENT
variables <- c("RC", "RT", "LC", "LT", "RKJC", "LKJC", "RTibiaToGlobal", "LTibiaToGlobal", "Force.Plate.1", "Force.Plate.2")

for (i in 1:length(DanielCSV)){
  DanielCSV[[i]] <- extractVariables(DanielCSV[[i]], variables)
}

for (i in 1:length(ThorCSV)){
  ThorCSV[[i]] <- extractVariables(ThorCSV[[i]], variables)
}
#-------------------------------------------------------------------------------------------------------------

#filter dataframe (mainly grw1 cause it doesn't seem to have been filtered in vicon)--------------------------
# for (i in 1:length(DanielCSV)){
#   DanielCSV[[i]] <- filterDataframe(DanielCSV[[i]])
# }
# 
# for (i in 1:length(ThorCSV)){
#   ThorCSV[[i]] <- filterDataframe(ThorCSV[[i]])
# }
#-------------------------------------------------------------------------------------------------------------

#80N threshold for Z force plate---------------------------------------
for (i in 1:length(DanielCSV)){
  DanielCSV[[i]] <- makeZero(DanielCSV[[i]])
}

for (i in 1:length(ThorCSV)){
  ThorCSV[[i]] <- makeZero(ThorCSV[[i]])
}
#----------------------------------------------------------------------

# #COP-------------------------------------------------------------------
for (i in 1:length(DanielCSV)){
  DanielCSV[[i]]$GRW1_mm_X <- COP(DanielCSV[[i]], leftorright = "left")$x
  DanielCSV[[i]]$GRW1_mm_Y <- COP(DanielCSV[[i]], leftorright = "left")$y
  DanielCSV[[i]]$GRW2_mm_X <- COP(DanielCSV[[i]], leftorright = "right")$x
  DanielCSV[[i]]$GRW2_mm_Y <- COP(DanielCSV[[i]], leftorright = "right")$y
}

for (i in 1:length(ThorCSV)){
  ThorCSV[[i]]$GRW1_mm_X <- COP(ThorCSV[[i]], leftorright = "left")$x
  ThorCSV[[i]]$GRW1_mm_Y <- COP(ThorCSV[[i]], leftorright = "left")$y
  ThorCSV[[i]]$GRW2_mm_X <- COP(ThorCSV[[i]], leftorright = "right")$x
  ThorCSV[[i]]$GRW2_mm_Y <- COP(ThorCSV[[i]], leftorright = "right")$y
}
----------------------------------------------------------------------

#Calculate FPA, SW-------------------------------------------------------------------------------------
for (i in 1:length(DanielCSV)){
  DanielCSV[[i]][is.na(DanielCSV[[i]])] <- 0
  DanielCSV[[i]]$FPA_Left <- FPA(DanielCSV[[i]])$Left
  DanielCSV[[i]]$FPA_Right <- FPA(DanielCSV[[i]])$Right
  DanielCSV[[i]]$StartLeftStep <- Step(DanielCSV[[i]])$StartLeft
  DanielCSV[[i]]$StartRightStep <- Step(DanielCSV[[i]])$StartRight
  DanielCSV[[i]]$EndLeftStep <- Step(DanielCSV[[i]])$EndLeft
  DanielCSV[[i]]$EndRightStep <- Step(DanielCSV[[i]])$EndRight
  DanielCSV[[i]]$SW <- StepWidth(DanielCSV[[i]])
  DanielCSV[[i]]$LposX <- PosVec(DanielCSV[[i]])$Lx
  DanielCSV[[i]]$LposY <- PosVec(DanielCSV[[i]])$Ly
  DanielCSV[[i]]$LposZ <- PosVec(DanielCSV[[i]])$Lz
  DanielCSV[[i]]$RposX <- PosVec(DanielCSV[[i]])$Rx
  DanielCSV[[i]]$RposY <- PosVec(DanielCSV[[i]])$Ry
  DanielCSV[[i]]$RposZ <- PosVec(DanielCSV[[i]])$Rz
  # DanielCSV[[i]]$LposX_V <- PosVec(DanielCSV[[i]])$Lx_V
  # DanielCSV[[i]]$LposY_V <- PosVec(DanielCSV[[i]])$Ly_V
  # DanielCSV[[i]]$LposZ_V <- PosVec(DanielCSV[[i]])$Lz_V
  # DanielCSV[[i]]$RposX_V <- PosVec(DanielCSV[[i]])$Rx_V
  # DanielCSV[[i]]$RposY_V <- PosVec(DanielCSV[[i]])$Ry_V
  # DanielCSV[[i]]$RposZ_V <- PosVec(DanielCSV[[i]])$Rz_V
  DanielCSV[[i]]$LKAM <- CalculateKam(DanielCSV[[i]], type = "actual")$Left
  DanielCSV[[i]]$RKAM <- CalculateKam(DanielCSV[[i]], type = "actual")$Right
  # DanielCSV[[i]]$LKAM_V <- CalculateKam(DanielCSV[[i]], type = "vicon")$Left
  # DanielCSV[[i]]$RKAM_V <- CalculateKam(DanielCSV[[i]], type = "vicon")$Right
  DanielCSV[[i]]$LKAM_norm <- (DanielCSV[[i]]$LKAM/1000)/(73*1.74)
  DanielCSV[[i]]$RKAM_norm <- (DanielCSV[[i]]$RKAM/1000)/(73*1.74)
}

for (i in 1:length(ThorCSV)){
  ThorCSV[[i]][is.na(ThorCSV[[i]])] <- 0
  ThorCSV[[i]]$FPA_Left <- FPA(ThorCSV[[i]])$Left
  ThorCSV[[i]]$FPA_Right <- FPA(ThorCSV[[i]])$Right
  ThorCSV[[i]]$StartLeftStep <- Step(ThorCSV[[i]])$StartLeft
  ThorCSV[[i]]$StartRightStep <- Step(ThorCSV[[i]])$StartRight
  ThorCSV[[i]]$EndLeftStep <- Step(ThorCSV[[i]])$EndLeft
  ThorCSV[[i]]$EndRightStep <- Step(ThorCSV[[i]])$EndRight
  ThorCSV[[i]]$SW <- StepWidth(ThorCSV[[i]])
  ThorCSV[[i]]$LposX <- PosVec(ThorCSV[[i]])$Lx
  ThorCSV[[i]]$LposY <- PosVec(ThorCSV[[i]])$Ly
  ThorCSV[[i]]$LposZ <- PosVec(ThorCSV[[i]])$Lz
  ThorCSV[[i]]$RposX <- PosVec(ThorCSV[[i]])$Rx
  ThorCSV[[i]]$RposY <- PosVec(ThorCSV[[i]])$Ry
  ThorCSV[[i]]$RposZ <- PosVec(ThorCSV[[i]])$Rz
  ThorCSV[[i]]$LKAM <- CalculateKam(ThorCSV[[i]])$Left
  ThorCSV[[i]]$RKAM <- CalculateKam(ThorCSV[[i]])$Right
  ThorCSV[[i]]$LKAM_norm <- (ThorCSV[[i]]$LKAM/1000)/(74*1.76)
  ThorCSV[[i]]$RKAM_norm <- (ThorCSV[[i]]$RKAM/1000)/(74*1.76)
}
#---------------------------------------------------------------------------------------------------

#Remove data from dataframe which is no longer needed-------------------------------------
dataframeColumns <- c('Force.Plate.1_Force.Z_N', 'Force.Plate.2_Force.Z_N', 'Force.Plate.1_COP.X_mm', 'Force.Plate.1_COP.Y_mm', 'Force.Plate.2_COP.X_mm', 'Force.Plate.2_COP.Y_mm', 'FPA_Left', 'FPA_Right', 'SW', 'StartLeftStep', 'StartRightStep', 'EndLeftStep', 'EndRightStep', 'LKAM_norm', 'RKAM_norm', 'LKAM', 'RKAM')

for (i in 1:length(DanielCSV)){
  DanielCSV[[i]] <- DanielCSV[[i]][,dataframeColumns]
}

for (i in 1:length(ThorCSV)){
  ThorCSV[[i]] <- ThorCSV[[i]][,dataframeColumns]
}
#-----------------------------------------------------------------------------------------

#Identify Peaks and save to list----------------------------------------------------------
DanielLeftPeak <- list()
DanielRightPeak <- list()
Peaks <- list()

for (i in 1:length(DanielCSV)){
  # Peaks <- PeakSelect(DanielCSV[[i]], leftorright = "left")
  # DanielLeftPeak[[i]] <- Peaks
  Peaks <- PeakSelect(DanielCSV[[i]], leftorright = "right")
  DanielRightPeak[[i]] <- Peaks
}

ThorLeftPeak <- list()
ThorRightPeak <- list()
Peaks <- list()

for (i in 1:length(ThorCSV)){
  Peaks <- PeakSelect(ThorCSV[[i]], leftorright = "left")
  ThorLeftPeak[[i]] <- Peaks
  Peaks <- PeakSelect(ThorCSV[[i]], leftorright = "right")
  ThorRightPeak[[i]] <- Peaks
}
# Plot <- plotSteps(DanielCSV[[1]], leftorright = "right")
# Plot
#-----------------------------------------------------------------------------------------

dir <- "C:\\Users\\bioengsu\\Desktop\\KAM\\Subject KAM CSV\\Trial 3\\"
setwd(dir)

DanielLeftFiles <- c()
DanielRightFiles <- c()
ThorLeftFiles <- c()
ThorRightFiles <- c()

for (i in 1:length(DanielFiles)){
  # DanielLeftFiles[i] <- paste("Daniel_Left", DanielFiles[i], sep = "_")
  DanielRightFiles[i] <- paste("Daniel_Right", DanielFiles[i], sep = "_")
}

for (i in 1:length(ThorFiles)){
  ThorLeftFiles[i] <- paste("Thor_Left", ThorFiles[i], sep = "_")
  ThorRightFiles[i] <- paste("Thor_Right", ThorFiles[i], sep = "_")
}

for (i in 1:length(DanielCSV)){
  # write.csv(peakDataframe(DanielCSV[[i]], DanielLeftPeak[[i]],leftorright = "left"), 
  #           file = DanielLeftFiles[i], row.names = FALSE, sep = ",")
  write.csv(peakDataframe(DanielCSV[[i]], DanielRightPeak[[i]],leftorright = "right"), 
            file = DanielRightFiles[i], row.names = FALSE, sep = ",")
}

for (i in 1:length(ThorCSV)){
  write.csv(peakDataframe(ThorCSV[[i]], ThorLeftPeak[[i]],leftorright = "left"), 
            file = ThorLeftFiles[i], row.names = FALSE, sep = ",")
  write.csv(peakDataframe(ThorCSV[[i]], ThorRightPeak[[i]],leftorright = "right"), 
            file = ThorRightFiles[i], row.names = FALSE, sep = ",")
}