#Cleaning second trial for Daniel and Thor 
#All inputs have been filtered in Vicon

rm(list = ls())

#functions---------------------------------------
determineIndices <- function(csvData){
  Indices <- c()
  
  for (i in 1:nrow(csvData)){
    if (csvData[i,1] == ""){
      Indices <- append(Indices, i)
    }
  }
  
  return(Indices)
}

#Store headers
header <- function(dataframe){
  head <- c()
  first <- ""
  second <- ""
  third <- ""
  for (i in 1:ncol(dataframe)){
    
    if (dataframe[1,i] != ""){
      first <- dataframe[1,i]
    }
    if (dataframe[2,i] != ""){
      second <- dataframe[2,i]
    }
    #third <- dataframe[3,i]
    
    if (first == ""){
      head <- append(head, paste(second, third, sep = "_"))
    } else {
      head <- append(head, paste(first, second, third, sep = "_"))
    }
    
  }
  return(head)
}

#Remove lines which aren't data
RemoveLines <- function(dataframe){
  dataframe <- dataframe[-(1:3),]
  
  return(dataframe)
}

#Downsample 10:1
Decimate <- function(dataframe){
  data2dec <- data.frame()
  
  #row index of first line of data
  index <- 4 #4
  
  rows <- seq(from=index, to=nrow(dataframe), by=10)
  
  data2dec <- dataframe[c(index-2,index-1, rows),-1]
  
  return(data2dec)
}

#------------------------------------------------

#read csv----------------------------------------
dir <- c("C:\\Users\\bioengsu\\Desktop\\KAM\\Subjects\\Daniel\\Trial 2", "C:\\Users\\bioengsu\\Desktop\\KAM\\Subjects\\Thor\\Trial 2")

setwd(dir[1])

DanielFiles <- list.files(dir[1])
ThorFiles <- list.files(dir[2])


DanielCSV <- lapply(DanielFiles, read.csv, header = FALSE, stringsAsFactors = FALSE, blank.lines.skip = FALSE, skip = 8)

setwd(dir[2])

ThorCSV <- lapply(ThorFiles, read.csv, header = FALSE, stringsAsFactors = FALSE, blank.lines.skip =  FALSE, skip = 8)
#------------------------------------------------

#Clean files-------------------------------------
#combine force with marker data

for (i in 1:length(DanielCSV)){
  Index <- determineIndices(DanielCSV[[i]])
  DanielCSV[[i]][1,] <- ""
  forceData <- DanielCSV[[i]][(Index[2]+8):(Index[6]-1), 1:25]

  DanielCSV[[i]] <- DanielCSV[[i]][1:(Index[2]-1),]
  
  forceData <- Decimate(forceData)

  DanielCSV[[i]] <- with(DanielCSV[[i]], cbind(DanielCSV[[i]], forceData))
  DanielCSV[[i]] <- DanielCSV[[i]][,-1]
}

for (i in 1:length(ThorCSV)){
  Index <- determineIndices(ThorCSV[[i]])
  ThorCSV[[i]][1,] <- ""
  forceData <- ThorCSV[[i]][(Index[2]+8):(Index[6]-1), 1:25]
  
  ThorCSV[[i]] <- ThorCSV[[i]][1:(Index[2]-1),]
  
  forceData <- Decimate(forceData)
  
  ThorCSV[[i]] <- with(ThorCSV[[i]], cbind(ThorCSV[[i]], forceData))
  ThorCSV[[i]] <- ThorCSV[[i]][,-1]
}
#------------------------------------------------

#headers for dataframes--------------------------
for (i in 1:length(DanielCSV)){
  colnames(DanielCSV[[i]]) <- header(DanielCSV[[i]])
  DanielCSV[[i]] <- RemoveLines(DanielCSV[[i]])
  
  for (j in 1:ncol(DanielCSV[[i]])){
    DanielCSV[[i]][,j] <- as.numeric(DanielCSV[[i]][,j])
  }
}

# for (j in 1:ncol(markerData)){
#   markerData[,j] <- as.numeric(markerData[,j])
# }

for (i in 1:length(ThorCSV)){
  colnames(ThorCSV[[i]]) <- header(ThorCSV[[i]])
  ThorCSV[[i]] <- RemoveLines(ThorCSV[[i]])
  
  for (j in 1:ncol(ThorCSV[[i]])){
    ThorCSV[[i]][,j] <- as.numeric(ThorCSV[[i]][,j])
  }
}
#------------------------------------------------

#write to csv------------------------------------
dir <- c("C:\\Users\\bioengsu\\Desktop\\KAM\\Cleaned CSV\\Trial 2\\Daniel", "C:\\Users\\bioengsu\\Desktop\\KAM\\Cleaned CSV\\Trial 2\\Thor")
setwd(dir[1])

for (i in 1:length(DanielFiles)){
  write.csv(DanielCSV[[i]], file = DanielFiles[i], row.names = FALSE)
}

setwd(dir[2])

for (i in 1:length(ThorFiles)){
  write.csv(ThorCSV[[i]], file = ThorFiles[i], row.names = FALSE)
}
#------------------------------------------------