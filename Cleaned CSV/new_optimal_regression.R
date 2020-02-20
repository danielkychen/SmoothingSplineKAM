#Script for creating new regression models of Trial 2 of Daniel and Thor

rm(list = ls())

library(car)
library(mgcv)

DanielFiles <- c(choose.files())
ThorFiles <- c(choose.files())

DanielCSV <- lapply(DanielFiles, read.csv)
ThorCSV <- lapply(ThorFiles, read.csv)

Daniel.df <- data.frame()
Thor.df <- data.frame()

#Combine CSV into one dataframe-------------------------------------
for (i in 1:length(DanielCSV)){
  Daniel.df <- rbind(Daniel.df, DanielCSV[[i]])
}

for (i in 1:length(ThorCSV)){
  Thor.df <- rbind(Thor.df, ThorCSV[[i]])
}
#-------------------------------------------------------------------

#Regression models--------------------------------------------------
stepRange <- seq(from=150, to=500, by=2)
fpaRange <- seq(from=-20, to=30, by = 1)

surfaceRegion <- data.frame(FPA = rep(fpaRange, each=length(stepRange)), StepWidth = rep(stepRange, times=length(fpaRange)))

Daniel_range.df <- Daniel.df[Daniel.df$FPA > -20,]

DanielScatter <- scatter3d(formula = Daniel_range.df$maxPeakKAM ~ Daniel_range.df$FPA + Daniel_range.df$StepWidth, 
                           surface = TRUE, fit = c("smooth"), xlab = "Foot Progression Angle", ylab = "KAM", zlab = "Step Width", 
                           axis.ticks = TRUE, model.summary = TRUE)

ThorScatter <- scatter3d(formula = Thor.df$maxPeakKAM ~ Thor.df$FPA + Thor.df$StepWidth, 
                         surface = TRUE, fit = c("smooth"), xlab = "Foot Progression Angle", ylab = "KAM", zlab = "Step Width", 
                         axis.ticks = TRUE, model.summary = TRUE)

Daniel_gaussMP <- gam(formula = Daniel_range.df$maxPeakKAM ~ s(Daniel_range.df$FPA, Daniel_range.df$StepWidth),
                    family = gaussian())

Daniel_gaussFP <- gam(formula = Daniel_range.df$firstPeakKAM ~ s(Daniel_range.df$FPA, Daniel_range.df$StepWidth),
                    family = gaussian())

Daniel_gaussSP <- gam(formula = Daniel_range.df$secondPeakKAM ~ s(Daniel_range.df$FPA, Daniel_range.df$StepWidth),
                    family = gaussian())

summary(Daniel_gaussMP)
summary(Daniel_gaussFP)
summary(Daniel_gaussSP)

ThorMP <- gam(formula = Thor.df$maxPeakKAM ~ s(Thor.df$FPA, Thor.df$StepWidth),
                      family = gaussian())

ThorFP <- gam(formula = Thor.df$firstPeakKAM ~ s(Thor.df$FPA, Thor.df$StepWidth),
                      family = gaussian())

ThorSP <- gam(formula = Thor.df$secondPeakKAM ~ s(Thor.df$FPA, Thor.df$StepWidth),
                      family = gaussian())

summary(ThorMP)
summary(ThorFP)
summary(ThorSP)
#-------------------------------------------------------------------