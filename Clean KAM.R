rm(list = ls())

library(signal)

#Script to clean up .csv files with marker and force plate data
#Functions-------------------------------------------------------------------

#Downsample 10:1
Decimate <- function(dataframe){
  data2dec <- data.frame()
  
  #row index of first line of data
  index <- match(1, dataframe[,1])
  
  rows <- seq(from=index, to=nrow(dataframe), by=10)
  
  data2dec <- dataframe[rows,-1]
  
  return(data2dec)
}

#Filter signal
filterSignal <- function(signal, cutoff, type, sampleFreq){
  
  bf <- butter(2, (2*cutoff)/(sampleFreq), type = type)
  y <- filtfilt(bf, signal)
  
  return(y)
}

#Store headers
header <- function(dataframe){
  head <- c()
  for (i in 1:ncol(dataframe)){
    
    if (dataframe[1,i] != ""){
      first <- dataframe[1,i]
      #print(first)
    }
    if (dataframe[2,i] != ""){
      second <- dataframe[2,i]
      #print(second)
    }
    if (dataframe[3,i] != ""){
      third <- dataframe[3,i]
      #print(third)
    }
    
    head <- append(head, paste(first, second, third, sep = "_"))
  }
  return(head)
}

#Remove lines which aren't data
RemoveLines <- function(dataframe){
  dataframe <- dataframe[-(1:3),]
  
  return(dataframe)
}


#----------------------------------------------------------------------------

directory <- "C:\\Users\\bioengsu\\Desktop\\KAM\\Subjects"

setwd(directory)

subjectDir <- list.dirs(directory)
subjectDir <- subjectDir[-1]

#List of subjects names------------------------------------------------------
subjects <- c()
for (i in 1:length(subjectDir)){
  string <- strsplit(subjectDir[i], split = '/')
  
  subjects <- append(subjects, string[[1]][2])
}
#----------------------------------------------------------------------------

#Read subject CSV file into data frame---------------------------------------
for (i in 1:length(subjects)){
  setwd(subjectDir[i])
  
  if (i == 1){
    AndrewBollen <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    AndrewBollenForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 2){
    AndrewWong <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    AndrewWongForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 3){
    Daniel <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    DanielForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 4){
    Frank <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    FrankForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 5){
    Jasmine <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    JasmineForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 6){
    Jonathan <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    JonathanForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 7){
    Kelly <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    KellyForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 8){
    Phoebe <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    PhoebeForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 9){
    Rakesh <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    RakeshForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 10){
    Thor <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    ThorForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } else if (i == 11){
    Wayne <- read.csv(list.files(getwd())[2], header = FALSE, stringsAsFactors = FALSE)
    WayneForce <- read.csv(list.files(getwd())[1], header = FALSE, stringsAsFactors = FALSE)
  } 
}

forceHeader <- AndrewBollenForce[1:3,2:ncol(AndrewBollenForce)]


#----------------------------------------------------------------------------

#Downsample Force 10:1-------------------------------------------------------
AndrewBollenForce <- Decimate(AndrewBollenForce)
AndrewWongForce <- Decimate(AndrewWongForce)
DanielForce <- Decimate(DanielForce)
FrankForce <- Decimate(FrankForce)
JasmineForce <- Decimate(JasmineForce)
JonathanForce <- Decimate(JonathanForce)
KellyForce <- Decimate(KellyForce)
PhoebeForce <- Decimate(PhoebeForce)
RakeshForce <- Decimate(RakeshForce)
ThorForce <- Decimate(ThorForce)
WayneForce <- Decimate(WayneForce)

AndrewBollenForce <- rbind(forceHeader, AndrewBollenForce)
AndrewWongForce <- rbind(forceHeader, AndrewWongForce)
DanielForce <- rbind(forceHeader, DanielForce)
FrankForce <- rbind(forceHeader, FrankForce)
JasmineForce <- rbind(forceHeader, JasmineForce)
JonathanForce <- rbind(forceHeader, JonathanForce)
KellyForce <- rbind(forceHeader, KellyForce)
PhoebeForce <- rbind(forceHeader, PhoebeForce)
RakeshForce <- rbind(forceHeader, RakeshForce)
ThorForce <- rbind(forceHeader, ThorForce)
WayneForce <- rbind(forceHeader, WayneForce)
#----------------------------------------------------------------------------

#Clean CSV-------------------------------------------------------------------
AndrewBollen <- with(AndrewBollen, AndrewBollen[,-1])
AndrewWong <- with(AndrewWong, AndrewWong[,-1])
Daniel <- with(Daniel, Daniel[,-1])
Frank <- with(Frank, Frank[,-1])
Jasmine <- with(Jasmine, Jasmine[,-1])
Jonathan <- with(Jonathan, Jonathan[,-1])
Kelly <- with(Kelly, Kelly[,-1])
Phoebe <- with(Phoebe, Phoebe[,-1])
Rakesh <- with(Rakesh, Rakesh[,-1])
Thor <- with(Thor, Thor[,-1])
Wayne <- with(Wayne, Wayne[,-1])
#----------------------------------------------------------------------------

#Headers---------------------------------------------------------------------
#String split to remove name
for (i in 1:ncol(AndrewBollen)){
  
  if (AndrewBollen[1,i] != ""){
    AndrewBollen[1,i] <- strsplit(AndrewBollen[1,i], split = ":")[[1]][2]
  }
  
  if (Wayne[1,i] != ""){
    Wayne[1,i] <- strsplit(Wayne[1,i], split = ":")[[1]][2]
  }
}

colnames(AndrewBollen) <- header(AndrewBollen)
colnames(AndrewWong) <- header(AndrewWong)
colnames(Daniel) <- header(Daniel)
colnames(Frank) <- header(Frank)
colnames(Jasmine) <- header(Jasmine)
colnames(Jonathan) <- header(Jonathan)
colnames(Kelly) <- header(Kelly)
colnames(Phoebe) <- header(Phoebe)
colnames(Rakesh) <- header(Rakesh)
colnames(Thor) <- header(Thor)
colnames(Wayne) <- header(Wayne)

colnames(AndrewBollenForce) <- header(AndrewBollenForce)
colnames(AndrewWongForce) <- header(AndrewWongForce)
colnames(DanielForce) <- header(DanielForce)
colnames(FrankForce) <- header(FrankForce)
colnames(JasmineForce) <- header(JasmineForce)
colnames(JonathanForce) <- header(JonathanForce)
colnames(KellyForce) <- header(KellyForce)
colnames(PhoebeForce) <- header(PhoebeForce)
colnames(RakeshForce) <- header(RakeshForce)
colnames(ThorForce) <- header(ThorForce)
colnames(WayneForce) <- header(WayneForce)

#Remove first few lines which were used for the headers
AndrewBollen <- RemoveLines(AndrewBollen)
AndrewWong <- RemoveLines(AndrewWong)
Daniel <- RemoveLines(Daniel)
Frank <- RemoveLines(Frank)
Jasmine <- RemoveLines(Jasmine)
Jonathan <- RemoveLines(Jonathan)
Kelly <- RemoveLines(Kelly)
Phoebe <- RemoveLines(Phoebe)
Rakesh <- RemoveLines(Rakesh)
Thor <- RemoveLines(Thor)
Wayne <- RemoveLines(Wayne)

AndrewBollenForce <- RemoveLines(AndrewBollenForce)
AndrewWongForce <- RemoveLines(AndrewWongForce)
DanielForce <- RemoveLines(DanielForce)
FrankForce <- RemoveLines(FrankForce)
JasmineForce <- RemoveLines(JasmineForce)
JonathanForce <- RemoveLines(JonathanForce)
KellyForce <- RemoveLines(KellyForce)
PhoebeForce <- RemoveLines(PhoebeForce)
RakeshForce <- RemoveLines(RakeshForce)
ThorForce <- RemoveLines(ThorForce)
WayneForce <- RemoveLines(WayneForce)
#----------------------------------------------------------------------------

#Combine Marker and Force Plate data-----------------------------------------
AndrewBollen <- cbind.data.frame(AndrewBollen, AndrewBollenForce)
AndrewWong <- cbind.data.frame(AndrewWong, AndrewWongForce)
Daniel <- cbind.data.frame(Daniel, DanielForce)
Frank <- cbind.data.frame(Frank, FrankForce)
Jasmine <- cbind.data.frame(Jasmine, JasmineForce)
Jonathan <- cbind.data.frame(Jonathan, JonathanForce)
Kelly <- cbind.data.frame(Kelly, KellyForce)
Phoebe <- cbind.data.frame(Phoebe, PhoebeForce)
Rakesh <- cbind.data.frame(Rakesh, RakeshForce)
Thor <- cbind.data.frame(Thor, ThorForce)
Wayne <- cbind.data.frame(Wayne, WayneForce)
#----------------------------------------------------------------------------

#Export data as CSV----------------------------------------------------------
setwd('C:\\Users\\bioengsu\\Desktop\\KAM\\Cleaned CSV')
MarkerList <- list(AndrewBollen, AndrewWong, Daniel, Frank, Jasmine, Jonathan, Kelly, Phoebe, Rakesh, Thor, Wayne)
MarkerNames <- c("Andrew Bollen.csv", "Andrew Wong.csv", "Daniel.csv", "Frank.csv", "Jasmine.csv", "Jonathan.csv",
                 "Kelly.csv", "Phoebe.csv", "Rakesh.csv", "Thor.csv", "Wayne.csv")

for (i in 1:length(MarkerList)){
  write.csv(MarkerList[[i]], file = MarkerNames[i], col.names = TRUE, row.names = FALSE)
}
#----------------------------------------------------------------------------



