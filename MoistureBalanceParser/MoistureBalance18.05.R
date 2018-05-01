cat("\n Moisture Balance data parser.  ", "Last edited: February 9, 2018.\n\n")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Parses data collected with Ohaus MB45 SN 1122272694
#
# Fit data to drying models taken from:
#  Sridhar, D. and Motappa Madhu, G. (2015) Drying Kinetics and Mathematical Modeling of Casuarina Equisetifolia
#     Wood Chips at Various Temperatures.  Period. Poly. Chem. Eng. 59(4) p288-295.
# 
# Data are fit using lm() or nls() methods with hard coded starting parameters for the nls() models
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# To execute type: source("MoistureBalance18.05.R", print.eval = TRUE)


# Load libraries

library(tcltk2)
library(ggplot2)
library(R2HTML)
library(tidyr)

# Clear all variables
# rm(list = ls())


# cat("\014")

setClass(Class="nlsfit",
         representation(
           fit="matrix",
           rsq="numeric",
           msg="character"
         )
)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Constants and initial values
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

modelResults <- NULL
allModelParams <- NULL
allModelFits <- NULL
allModelFitsLong <- NULL
parserID = "MB"
parserVersion = "18.05"
parserCounter <- 0
sumData <- NULL
logTime <- Sys.time()
todaysdate <- date()
modelAbbr <- c("Ex","HP", "Le", "Li", "Lo", "Loga", "MHP", "MP", "MPII", "Pa", "SF", "Tt", "WS")
modelParameterList <- c("a", "b", "c", "d", "e", "k", "ka", "kb", "kc", "kd", "ke", "n", "m", "p", "q", "D", "L", "V", "T")
modelEstimates <- data.frame(matrix(data = NA, nrow = 1, ncol = length(modelParameterList)))
colnames(modelEstimates) <- modelParameterList
modelPr.gt.ts <- modelEstimates
modelStdErrors <- modelEstimates
modeltValues <- modelEstimates
nModels = 13  # Plot and report model results from the top nModels based on cor(moisture, fit)
nModelResults = 13 

modelResultsList <- c("date", "time", "datafile", "sample", "dryingProfile", "finalTemp", "moisture", "modelAbbr", "modelName", "correlation",
                      modelParameterList,
                      paste("StdErr_", modelParameterList, sep = ""),
                      paste("tValue_", modelParameterList, sep = ""),
                      paste("Pr.gt.t_", modelParameterList, sep = ""),
                      "parserVersion", "fitError"
)

modelNames <- c("Exponential Fit", "Henderson Pabis Fit", "Lewis Fit", "Linear Fit", "Logarithmic Fit", "Logarithmic Fit (Yagcioglu et al. 1999)",
                "Modified Henderson Pabis Fit", "Modified Page Fit", "Modified Page II Fit", "Page Fit", "Simplified Fick Diffusion Fit",
                "Two-Term Fit (Henderson 1974)", "Wang and Singh Fit")
modelFormula <- c("log(moisture) ~ a * time + b", "moisture ~ a * exp(-k * time)", "moisture ~ exp(-k * time)", "moisture ~ a * time + b", "exp(moisture) ~ a * time + b",
                  "moisture ~ a * exp(-k * time) + c", "moisture ~ a * exp(-ka * time) + b * exp(-kb * time) + c * exp(-kc * time)", "moisture ~ exp(-(k * time) ^ n)",
                  "moisture ~ exp(-k * (time / L ^ 2) ^ n)", "moisture ~ exp(-k * time ^ n)", "moisture ~ a * exp((-c * time) / L ^ 2)",
                  "moisture ~ a * exp(-ka * time) + b * exp(-kb * time)", "moisture ~ 1 + a * time + b * time ^ 2")
radiobuttondone <- tclVar(0)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function chooseMultiple allows user to choose between single data file and
# whole direcory analysis
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

chooseMultiple = function(){
  rbVal = as.numeric(tclvalue(rbValue))
  tkdestroy(tt)
  if (rbVal==1)
    tkmessageBox(title = "Select folder", message="Analyze all files in selected directory.")
  if (rbVal==2)
    tkmessageBox(title = "Select file(s)", message="Analyze selected file(s) only.")
  tclvalue(radiobuttondone) <- 1
  ms <<- rbVal
}

decimalTime = function(time) {
  time = as.numeric(unlist(strsplit(time, ":")))
  time = time[1]*60+time[2]+time[3]/60
  return(time)
}

fitData <- function() {
  
  # Generate data to test fitting models
  
  # Ex <-  Exp_fit(rawData$dectime, rawData$moisture) 
  
  t <- seq(0, 10, 0.2)
  m <- exp(-0.4 * t)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # HP <-  HendersonPabis_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0, 10, 0.2)
  m <- 0.8 * exp(-0.5 * t)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # Le <-  Lewis_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0, 10, 0.2)
  m <- exp(-0.4 * t)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # Li <-  Linear_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0, 10, 0.2)
  m <- -0.08 * t + 0.88
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # Lo <- Log_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0.2, 10, 0.2)
  m <- log(5 * t)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # Loga <- Logarithmic_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0, 10, 0.2)
  m <- 0.4 * exp(-0.5 * t) + 0.4
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # MHP <- ModifiedHendersonPabis_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0.2, 10, 0.2)
  m <- 0.3 * exp(-0.2 * t) + 0.2 * exp(-0.5 * t) + 0.4 * exp(-0.3 * t)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # MP <- ModifiedPage_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0, 10, 0.2)
  m <- exp(-(0.5 * t) ^ 0.5)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # MPII <- ModifiedPageII_fit(rawData$dectime, rawData$moisture, 0.1)
  
  k = 0.14
  L = 0.43
  n = 1.29
  t <- seq(0, 10, 0.2)
  m <- exp(-(k * (t / L ^ 2) ^ n))
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # Pa <- Page_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0, 10, 0.2)
  m <- exp(-0.5 * t ^ 0.6)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # SF <- SimplifiedFick_fit(rawData$dectime, rawData$moisture, 1)
  
  L = 0.6
  t <- seq(0, 10, 0.2)
  m <- 0.8 * exp((-0.5 * time) / L ^ 2)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # Tt <- TwoTerm_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0, 10, 0.2)
  m <- 0.3 * exp(-0.5 * time) + 0.5 * exp(-0.1 * time)
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
  
  # WS <- WangSingh_fit(rawData$dectime, rawData$moisture)
  
  t <- seq(0, 10, 0.2)
  m <-  1 - 0.5 * time + 0.05 * time ^ 2
  m <- m + rnorm(length(t), 0, 0.1 * sd(m))
  plot(t, m)
  time <- t
  moisture <- m
}

# ModifiedPageII_fit(time, moisture)

Exp_fit <- function(time, moisture) {
  nonZero <- which(moisture > 0)
  moistureNZ <- moisture[nonZero]
  time <- time[nonZero]
  tryResult <- tryCatch(lm(log(moistureNZ) ~ time),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "lm") {
    model <- lm(log(moistureNZ) ~ time)
    Ex <- exp(predict(model))
    # plot(timeNZ, moistureNZ, pch = 16, cex = 1.3, col = "blue", main = "Exponential Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(timeNZ, Ex, lty=2, col="red", lwd=3)
    t <- time
    w <- Ex
    rsq <- as.numeric(cor(moistureNZ, Ex))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

HendersonPabis_fit <- function(time, moisture) {
  tryResult <- tryCatch(nls(moisture ~ a * exp(-k * time), start = list(a = 0.8, k = 0.5)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {
    model <- nls(moisture ~ a * exp(-k * time), start = list(a = 0.8, k = 0.5))
    HP <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Henderson Pabis Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, HP,lty=2,col="red",lwd=3)
    t <- time
    w <- HP
    rsq <- as.numeric(cor(moisture, HP))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

Lewis_fit <- function(time, moisture) {
  tryResult <- tryCatch(nls(moisture ~ exp(-k * time), start = list(k = 0.4)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {
    model <- nls(moisture ~ exp(-k * time), start = list(k = 0.4))
    Le <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Lewis Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, Le, lty=2, col="red", lwd=3)
    t <- time
    w <- Le
    rsq <- as.numeric(cor(moisture, Le))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

Linear_fit <- function(time,moisture) {
  tryResult <- tryCatch(lm(moisture ~ time),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "lm") {
    model <- lm(moisture ~ time)
    Li <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Linear Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, Li, lty=2, col="red", lwd=3)
    t <- time
    w <- Li
    rsq <- as.numeric(cor(moisture, Li))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

Log_fit <- function(time,moisture) {
  tryResult <- tryCatch(lm(exp(moisture) ~ time),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "lm") {
    model <- lm(exp(moisture) ~ time)
    #modelFinite <- which(predict(model) > 0)
    #Lo <- log(predict(model)[modelFinite])
    Lo <- log(predict(model))
    # plot(time, moisture, pch = 16, cex = 1.3, col = "blue", main = "Logarithmic Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, Lo, lty=2, col="red", lwd=3)
    t <- time
    w <- Lo
    rsq <- as.numeric(cor(moisture, Lo))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

Logarithmic_fit <- function(time, moisture) {
  tryResult <- tryCatch(nls(moisture ~ a * exp(-k * time) + c, start = list(a = 0.4, k = 0.5, c = 0.4)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {
    model <- nls(moisture ~ a * exp(-k * time) + c, start = list(a = 0.4, k = 0.5, c = 0.4))
    Loga <- predict(model) 
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Logarithmic Fit (Yagcioglu et al. 1999)", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, Loga, lty=2, col="red", lwd=3)
    t <- time
    w <- Loga
    rsq <- as.numeric(cor(moisture, Loga))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

ModifiedHendersonPabis_fit <- function(time,moisture) {
  
  nonZero <- which(moisture > 0)
  moistureNZ <- moisture[nonZero]
  time <- time[nonZero]
  
  tryResult <- tryCatch(nls(moistureNZ ~ a * exp(-ka * time) + b * exp(-kb * time) + c * exp(-kc * time), start = list(a = 0.3, b = 0.2, c = 0.4, ka = 0.2, kb = 0.5, kc = 0.3)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {
    model <- nls(moistureNZ ~ a * exp(-ka * time) + b * exp(-kb * time) + c * exp(-kc * time), start = list(a = 0.3, b = 0.2, c = 0.4, ka = 0.2, kb = 0.5, kc = 0.3))
    MHP <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Modified Henderson Pabis Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, MHP,lty=2,col="red",lwd=3)
    t <- time
    w <- MHP
    rsq <- as.numeric(cor(moistureNZ, MHP))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

ModifiedPage_fit <- function(time, moisture) {
  tryResult <- tryCatch(nls(moisture ~ exp(-(k * time) ^ n), start = list(k = 0.5, n = 0.5)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {
    model <- nls(moisture ~ exp(-(k * time) ^ n), start = list(k = 0.5, n = 0.5))
    MP <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Modified Page Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, MP,lty=2,col="red",lwd=3)
    t <- time
    w <- MP
    rsq <- as.numeric(cor(moisture, MP))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

ModifiedPageII_fit <- function(time, moisture) {
  tryResult <- tryCatch(nls(moisture ~ exp(-(k * (time / L ^ 2) ^ n)), start = list(k = 0.5, L = 0.1, n = 0.5)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {
    model <- nls(moisture ~ exp(-(k * (time / L ^ 2) ^ n)), start = list(k = 0.5, L = 0.1, n = 0.5))
    MPII <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = paste("Modified Page II Fit (L = ", L, ")", sep = ""), xlab = "Time (min)", ylab = "Moisture")
    # lines(time, MPII, lty=2, col="red", lwd=3)
    t <- time
    w <- MPII
    rsq <- as.numeric(cor(moisture, MPII))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

Page_fit <- function(time,moisture) {
  tryResult <- tryCatch(nls(moisture ~ exp(-k * time ^ n), start = list(k = -0.5, n = 0.6)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {
    model <- nls(moisture ~ exp(-k * time ^ n), start = list(k = -0.5, n = 0.6))
    Pa <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Page Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, Pa,lty=2,col="red",lwd=3)
    t <- time
    w <- Pa
    rsq <- as.numeric(cor(moisture, Pa))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

SimplifiedFick_fit <- function(time, moisture) {
  tryResult <- tryCatch(nls(moisture ~ a * exp((-c * time) / L ^ 2), start = list(a = 0.8, c = 0.5, L = 1.0)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {
    model <- nls(moisture ~ a * exp((-c * time) / L ^ 2), start = list(a = 0.8, c = 0.5, L = 1.0))
    SF <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = paste("Simplified Fick Diffusion Fit (L = ", L, ")", sep = ""), xlab = "Time (min)", ylab = "Moisture")
    # lines(time, SF,lty=2,col="red",lwd=3)
    t <- time
    w <- SF
    rsq <- as.numeric(cor(moisture, SF))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

TwoTerm_fit <- function(time, moisture) {
  tryResult <- tryCatch(nls(moisture ~ a * exp(-ka * time) + b * exp(-kb * time), start = list(a = 0.3, ka = 0.5, b = 0.5, kb = 0.1)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {  
    model <- nls(moisture ~ a * exp(-ka * time) + b * exp(-kb * time), start = list(a = 0.3, ka = 0.5, b = 0.5, kb = 0.1))
    Tt <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Two Term Fit (Henderson 1974)", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, Tt, lty=2, col="red", lwd=3)
    t <- time
    w <- Tt
    rsq <- as.numeric(cor(moisture, Tt))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

WangSingh_fit <- function(time, moisture) {
  tryResult <- tryCatch(nls(moisture ~ 1 + a * time + b * time ^ 2, start = list(a = -0.5, b = 0.05)),
                        error=function(e) {
                          print(e)
                          print(strsplit(as.character(e), ":")[[1]][2])
                        }
  )
  if (class(tryResult) == "nls") {  
    model <- nls(moisture ~ 1 + a * time + b * time ^ 2, start = list(a = -0.5, b = 0.05))
    WS <- predict(model)
    # plot(time,moisture, pch = 16, cex = 1.3, col = "blue", main = "Wang and Singh Fit", xlab = "Time (min)", ylab = "Moisture")
    # lines(time, WS, lty=2, col="red", lwd=3)
    t <- time
    w <- WS
    rsq <- as.numeric(cor(moisture, WS))
    err <- NA
  } else {
    model <- NA
    t <- NA
    w <- NA
    rsq <- NA
    err <- gsub("[\r\n]", "", tryResult)
  }
  return(list(model, t, w, rsq, err))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Select folder or single file to be processed
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

ms = NA
tt <- tktoplevel()
tktitle(tt) <- "Select folder of selected file(s)"
rb1 <- tkradiobutton(tt)
rb2 <- tkradiobutton(tt)
rbValue <- tclVar(1)
tkconfigure(rb1,variable=rbValue,value=1)
tkconfigure(rb2,variable=rbValue,value=2)
tkgrid(tklabel(tt,text="Select folder or selected file(s) analysis:"))
tkgrid(tklabel(tt,text="Select folder"),rb1)
tkgrid(tklabel(tt,text="Select file(s)"),rb2)
OK.but <- tkbutton(tt, text="OK", command = chooseMultiple)
tkgrid(OK.but)
tkfocus(tt)
tkwait.variable(radiobuttondone)

if (ms == 1) {
  
  # Analyze folder
  
  dataFolder <- choose.dir(default = "C:/", caption = "Select folder")
  dataFolder <- gsub("\\\\", "/", dataFolder)
  dirList <- dir(dataFolder, recursive = TRUE)
  
  # Select .txt files only
  
  xxx.txt <- dirList[grep(".txt", dirList)]
} else {
  
  # Analyze selected individual file or files
  
  infile = choose.files(paste("C:/", "*.txt", sep = ""), filters = c("txt","All"),
                        caption = "Selecte file(s) to be analyzed.")
  infile <- gsub("\\\\", "/", infile)
  fName <- paste("/", strsplit(infile, "/")[[1]][length(strsplit(infile, "/")[[1]])], sep = "")
  dataFolder <- strsplit(infile, fName)[[1]][1]
  xxx.txt <- sapply(strsplit(infile, dataFolder), function(x) x[[2]])
  xxx.txt <- gsub("/", "", xxx.txt)
}

# Ignore the _info.txt files

infoCount <- grep("_info.txt", xxx.txt)
if (length(infoCount) > 0 ) {
  xxx.txt <- xxx.txt[-infoCount]
}



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Select folder containing data to be parsed.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

logTime <- gsub(":", "", logTime)
# dataFolder <- choose.dir(default = "C:/", caption = "Select folder")
# dataFolder <- gsub("\\\\", "/", dataFolder)

logfName <- paste(strsplit(dataFolder, "/")[[1]][length(strsplit(dataFolder, "/")[[1]])], "_", parserID, parserVersion, "_", logTime, ".log", sep = "")
logfName <- paste(dataFolder, logfName, sep = "/")
infoName <-  paste(strsplit(dataFolder, "/")[[1]][length(strsplit(dataFolder, "/")[[1]])], parserID, parserVersion, "_", logTime, "_info.txt", sep = "")
infoName <- paste(dataFolder, infoName, sep = "/")

# Get computer name

nodename <- Sys.info()["nodename"]
cat("Computer name:", nodename, "\n", file = logfName, append = TRUE)
cat("Parser version:", parserVersion, "\n", file = logfName, append = TRUE)
cat("Sequence name:", dataFolder, "\n", file = logfName, append = TRUE)
cat("Date:", todaysdate, "\n", file = logfName, append = TRUE)
cat("\n\nSample info written to:", infoName, "\n\n", file = logfName, append = TRUE)
cat("\n\nMoistureBalance parsed data written to:\n\n", file = logfName, append = TRUE)

nodename <- Sys.info()["nodename"]
cat("Computer name:", nodename, "\n", file = infoName, append = TRUE)
cat("Parser version:", parserVersion, "\n", file = infoName, append = TRUE)
cat("Sequence name:", dataFolder, "\n", file = infoName, append = TRUE)
cat("Date:", todaysdate, "\n", file = infoName, append = TRUE)
cat("\n\nParser log data written to :", logfName, "\n\n", file = infoName, append = TRUE)
cat("infile, testID, switchOffMode, dryingProfile, dryingTemp, resultUnits, initialWeight, finalWeight, elapsedTime, finalPanTemp, finalResult, sampleNumber\n", file = infoName, append = TRUE) 



for (f in xxx.txt) {
  
  # Get file name
  
  fName <- strsplit(f, "/")[[1]][length(strsplit(f, "/")[[1]])]
  csvName <- sub(".txt", paste("_", parserID, parserVersion,".csv", sep = ""), f)
  modName <- sub(".txt", paste("_", parserID, parserVersion,".mod", sep = ""), f)
  infile <- paste(dataFolder, f, sep = "/")
  if (file.info(infile)$size > 0) {
    MBData <- read.csv(infile, sep = "\t", header = FALSE, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
    
    # There may still be other .txt files in the directory that do not contain Ohaus MB45 data.
    # Look for specific text in line 2 of the .txt file.  Hopefully this elimiates non-data files.
    
    if (MBData[2,] == "OHAUS MB45 SN   1122272694") {
      testID <- grep("TEST ID", MBData[,], ignore.case = TRUE)  
      switchOffMode <- grep("Switchoff Mode", MBData[,], ignore.case = TRUE)  
      dryingProfile <- grep("Drying Profile", MBData[,], ignore.case = TRUE)  
      dryingTemp <- grep("Drying Temp", MBData[,], ignore.case = TRUE)  
      resultUnits <- grep("Result Units", MBData[,], ignore.case = TRUE)  
      initialWeight <- grep("Initial Weight", MBData[,], ignore.case = TRUE)
      finalWeight <- grep("Final Weight", MBData[,], ignore.case = TRUE)
      elapsedTime <- grep("Elapsed Time", MBData[,], ignore.case = TRUE)  
      finalPanTemp <- grep("Final Pan Temp", MBData[,], ignore.case = TRUE)  
      finalResult <- grep("Final Result", MBData[,], ignore.case = TRUE)
      endFlag <- grep("--End--", MBData[,], ignore.case = TRUE)  
      nSamples <- length(endFlag)
      
      # Verify that data file is correctly formatted  
      
      if (length(initialWeight) == 2 * nSamples) {
        
        # Remove duplicate Initial Weight values
        
        initialWeight <- initialWeight[seq(1, by = 2, len = nSamples)]
        beginData <- initialWeight + 1
        endData <- elapsedTime - 1
        
        for (i in 1:nSamples) {
          
          # Verify that there are more than 2 data points in the data file
          
          if (endData[i] - beginData[i] > 1) {
            
            # Proceed with data processing if there are 2 or more data points
            
            
            predictDF <- NULL
            predictDFW <- NULL
            predictNext <- NULL
            corr <- NULL
            
            #   Extract parameter info
            
            testID[i] <- gsub("\\s+", "", strsplit(MBData[testID[i],], ":")[[1]][2])
            switchOffMode[i] <- gsub("\\s+", "", strsplit(MBData[switchOffMode[i],], "Mode.")[[1]][2])
            dryingProfile[i] <- gsub("\\s+", "", strsplit(MBData[dryingProfile[i],], "Profile.")[[1]][2])  
            dryingTemp[i] <- gsub("\\s+", "", strsplit(MBData[dryingTemp[i],], "Temp.")[[1]][2])
            resultUnits[i] <- gsub("\\s+", "", strsplit(MBData[resultUnits[i],], "Units.")[[1]][2])
            
            initialWeight[i] <- gsub("\\s+", "", strsplit(MBData[initialWeight[i],], "Weight.")[[1]][2])
            elapsedTime[i] <- gsub("\\s+", "", strsplit(MBData[elapsedTime[i],], "Time.")[[1]][2])
            finalWeight[i] <- gsub("\\s+", "", strsplit(MBData[finalWeight[i],], "Weight.")[[1]][2])
            
            initialWeight[i] <- gsub("[^0-9\\.]", "", initialWeight[i])
            finalWeight[i] <- gsub("[^0-9\\.]", "", finalWeight[i])
            finalPanTemp[i] <- gsub("\\s+", "", strsplit(MBData[finalPanTemp[i],], "Temp.")[[1]][2])
            finalPanTemp[i] <- gsub("[^0-9\\.]", "", finalPanTemp[i])
            finalResult[i] <- gsub("\\s+", "", strsplit(MBData[finalResult[i],], "Result.")[[1]][2])
            finalResult[i] <- gsub("[^0-9\\.]", "", finalResult[i])
            
            rawData <- MBData[beginData[i]:endData[i],]
            rawData <- gsub("[^0-9\\.\\:]", " ", rawData)
            rawData <- gsub("\\s+", ",", rawData)
            rawData <- paste(rawData, i)
            rawData <- read.table(text = rawData, sep = ",", header = FALSE,
                                  col.names = c('time', 'temperature', 'moisture', 'sample'),
                                  as.is = TRUE)
            
            finalTemp <- as.numeric(finalPanTemp[i])
            moisture <- as.numeric(finalResult[i])
            
            
            # Convert time from character string to time element.  This will append the current date to the date-time element.
            # The time may need to be converted to decimal time depending on how the fitting method represents time.
            
            # rawData$time <- as.POSIXct(rawData$time,format="%H:%M:%S")
            
            dectime = sapply(rawData$time,decimalTime)
            
            # Calculate moisture ratio MR, and include in rawData
            
            MR <- round(1 - rawData$moisture / max(rawData$moisture, na.rm = TRUE), 3)
            
            # Calculate masses based on initialWeight(s)
            
            mass <- as.numeric(initialWeight[i]) - as.numeric(initialWeight[i]) * rawData$moisture / 100
            
            # Include mass in rawData
            
            rawData <- cbind(rawData[, 1:3], dectime, mass, MR, rawData[, 4])
            
            # For some reason this process removes the name from the sample columm.  Re-assign column names.
            
            names(rawData) <-  c('time', 'temperature', 'moisture', 'dectime', 'mass', 'MR', 'sample')
            
            # Fit data to models
            
            Ex <-  Exp_fit(rawData$dectime, rawData$MR) 
            HP <-  HendersonPabis_fit(rawData$dectime, rawData$MR)
            Le <-  Lewis_fit(rawData$dectime, rawData$MR)
            Li <-  Linear_fit(rawData$dectime, rawData$MR)
            Lo <- Log_fit(rawData$dectime, rawData$MR)
            Loga <- Logarithmic_fit(rawData$dectime, rawData$MR)     
            MHP <- ModifiedHendersonPabis_fit(rawData$dectime, rawData$MR)           
            MP <- ModifiedPage_fit(rawData$dectime, rawData$MR) 
            MPII <- ModifiedPageII_fit(rawData$dectime, rawData$MR)          
            Pa <- Page_fit(rawData$dectime, rawData$MR) 
            SF <- SimplifiedFick_fit(rawData$dectime, rawData$MR)
            Tt <- TwoTerm_fit(rawData$dectime, rawData$MR)
            WS <- WangSingh_fit(rawData$dectime, rawData$MR)
            
            # Tabulate time and predicted data for ALL models. Fill "missing" rows with NAs. Append raw data.
            
            # Use get() to access variable with variable name string.
            
            # Generate wide data frame of predicted data. Fill "missing data with NAs. Append raw data.
            # Previous versions of this code used the sorted corrDF to generate predictDFW so that the
            # models were added in decending order of correlation.  In order to append the model data to
            # the rawData for all data files, the model data is now added in order according to modelAbbr[].
            
            # Initialize predictDFW[] with time and MR from rawData[]
            
            predictDFW <- data.frame(t = rawData$dectime, raw = rawData$MR)            
            
            # Check for convergent models based on Pearson correlation coefficients.  If the model did not converge, the correlation is set to NA.
            
            for (j in 1:length(modelAbbr)) {
              corr <- c(corr, round(as.numeric((get(modelAbbr[j])[[4]])), 8))
            }
            corrDF <- data.frame(corr, modelAbbr, modelNames, stringsAsFactors = FALSE)
            nModels <- sum(!is.na(corrDF$corr))
            if (nModels > 0) {
              
              # Proceed with data processing to generate predictDFW[]
              
              # Sort corrDF by Pearson correlation, with NAs last.
              
              corrDF <- corrDF[order(corrDF$corr, decreasing = TRUE),]
              
              for (k in 1:length(modelAbbr)) {
                if (!is.na(get(modelAbbr[k])[[2]][1])) {
                  t <- get(modelAbbr[k])[[2]]
                  
                  # As of this version the model output contains the predicted moisture values in position 3 of the list,
                  # so it is no longer necessary to calculate using predict().
                  
                  # p <- predict(get(modelAbbr[k])[[1]])
                  
                  pred <- get(modelAbbr[k])[[3]]
                  predictNext <- data.frame(t, pred)
                  colnames(predictNext)[2] <- modelAbbr[k]
                  rowPredDFW <- nrow(predictDFW)
                  rowPredNext <- nrow(predictNext)
                  if (rowPredDFW == rowPredNext) {
                    
                    # rowPredNext is ready to add to predictDFW without alterations.
                    
                    #predictDFW <- cbind(predictDFW, predictNext[2])
                  } else {
                    if (rowPredNext > rowPredDFW) {
                      m <- data.frame(matrix(NA, nrow = (rowPredNext - rowPredDFW), ncol = ncol(predictDFW)))
                      names(m) <- names(predictDFW)
                      predictDFW <- rbind(predictDFW, m)
                    } else {
                      
                      # rowPredDFW > rowPredNext
                      
                      m <- data.frame(matrix(NA, nrow = (rowPredDFW - rowPredNext), ncol = ncol(predictNext)))
                      names(m) <- names(predictNext)
                      predictNext <- rbind(predictNext, m)
                    }
                  }
                } else {
                  
                  # No model data for modelAbbr[k].
                  
                  predictNext <- data.frame(matrix(NA, nrow = rowPredDFW, ncol = 2))
                  colnames(predictNext)[2] <- modelAbbr[k]
                }
                predictDFW <- cbind(predictDFW, predictNext[2])
              }
              
            } else {
              
              # No models fit the data.
              
              predictNext <- data.frame(matrix(NA, nrow = nrow(predictDFW), ncol = length(modelAbbr)))
              colnames(predictNext) <- modelAbbr
              predictDFW <- cbind(predictDFW, predictNext)
            }
            
            # Tabulate long data frame with data for top nModels for plotting.
            
            predictDF <- gather(predictDFW, model, MR, -t)
            tMax <- max(predictDF$t, na.rm = TRUE)
            mMax <- max(predictDF$MR, na.rm = TRUE)          
            
            for (q in 1:nModels) {
              if (q == 1) {
                lFinal <- paste("\n\ncorr_", corrDF$modelAbbr[q], " = ", round(as.numeric(corrDF$corr[q]), 4), sep = "")
              } else {
                lFinal <- paste(lFinal, "\n ", "corr_", corrDF$modelAbbr[q], " = ", round(as.numeric(corrDF$corr[q]), 4), sep = "")
              }
            }  #  Not sure the for loop in q needs to end here or not.
            
            plotTitle <- paste("Best fit for ", fName, if(nSamples > 1) {paste("[", i, "]", sep = "")}, sep = "")
            
            # Set .png file name.
            
            if (nSamples > 1) {
              pngName <- sub(".txt", paste("_", i, "_", parserID, parserVersion,".png", sep = ""), f)
              modelcsv <- sub(".txt", paste("_", i, "_", parserID, parserVersion,"_model.csv", sep = ""), f)
            } else {
              pngName <- sub(".txt", paste("_", parserID, parserVersion,".png", sep = ""), f)
              modelcsv <- sub(".txt", paste("_", parserID, parserVersion,"_model.csv", sep = ""), f)
            }
            
            HTMLoutput <- paste(dataFolder, sub(".csv", "Rep.html", modelcsv), sep = "/")
            
            # Plot raw data and top models
            
            plotData <- data.frame(t = rawData$dectime, MR = rawData$MR, model = "raw")
            png(paste(dataFolder, "/", pngName, sep = ""))
            g <- ggplot(data=plotData, aes(x=t, y=MR)) +
              geom_point(size = 4) +
              geom_line(data=predictDF, aes(x=t, y=MR, color = model), size = 1.5) +
              ggtitle(plotTitle) +
              annotate("text", label = lFinal, x = 0.9 * tMax, y = 1.2 * mMax, size = 5, hjust = "right") 
            print(g)
            dev.off()
            
            HTML("", append = FALSE, file = HTMLoutput)
            HTML.title(paste("Regression summary info for ", fName, sep = ""), Align = "center", HR = 3, file = HTMLoutput)
            HTML.title(Sys.time(), Align = "center", HR = 4, file = HTMLoutput)
            HTML.title(paste("Parser version:", parserVersion), Align = "center", HR = 4, file = HTMLoutput)
            HTML("<hr>",file = HTMLoutput)
            HTMLInsertGraph(paste(dataFolder, "/", pngName, sep = ""), file = HTMLoutput, GraphBorder = 3, Align = "center")
            
            for (p in 1:nModelResults) {
              HTML(paste("Model summary for the ",  corrDF$modelNames[p]), file = HTMLoutput)
              print(paste("Model summary for the ",  corrDF$modelNames[p]))
              
              # For some reason this stopped working...
              HTML(summary(get(corrDF$modelAbbr[p])[[1]]), file = HTMLoutput)
              
              # Extract and collate parameters from nls() and lm() fits
              
              print(paste("Model data for", corrDF$modelAbbr[p]))
              
              modelEquation <- modelFormula[corrDF$modelAbbr[p] == modelAbbr]
              modelEstimates[] <- NA
              modelStdErrors <- modelEstimates
              modeltValues <- modelEstimates
              modelPr.gt.ts <- modelEstimates
              
              
              if (is.na(get(corrDF$modelAbbr[p])[[1]][1])) {
                
                # nls() reported model fitting error. Could be a result of zero moisture. Otherwise results left in initial state of NA.
                # Error message reported.
                
              } else {
                
                modelParams <- summary(get(corrDF$modelAbbr[p])[[1]])$coefficients
                switch(corrDF$modelAbbr[p][[1]],
                       
                       # models generated by lm()
                       
                       "Ex" = {
                         print("Exponential fit found...")
                         rownames(modelParams) <- c("b", "a")
                       },
                       
                       "Li" = {
                         print("Linear fit found...")
                         rownames(modelParams) <- c("b", "a")
                       },
                       
                       "Lo" = {
                         print("Log function found...")
                         rownames(modelParams) <- c("b", "a")
                       },
                       
                       {
                         
                         # models generated by nls()
                         
                       }
                )
                modelParamsInfo <- colnames(modelParams)
                modelParamsNames <- rownames(modelParams)
                modelEstimates[modelParamsNames] <- modelParams[modelParamsNames, modelParamsInfo[1]]
                modelStdErrors[modelParamsNames] <- modelParams[modelParamsNames, modelParamsInfo[2]]
                modeltValues[modelParamsNames] <- modelParams[modelParamsNames, modelParamsInfo[3]]
                modelPr.gt.ts[modelParamsNames] <- modelParams[modelParamsNames, modelParamsInfo[4]]
              }
              
              colnames(modelStdErrors) <- paste("StdError_", colnames(modelStdErrors), sep = "")
              colnames(modeltValues) <- paste("tValue_", colnames(modeltValues), sep = "")
              colnames(modelPr.gt.ts) <- paste("Pr.gt.t_", colnames(modelPr.gt.ts), sep = "")
              
              # This generates a dataframe with rownames = model parameter names
              
              sampleID <- i
              dryingTime <- elapsedTime[i]
              initialTemp <- rawData$temperature[1]
              finalTemp <- as.numeric(finalPanTemp[i])
              sampleMass <- as.numeric(finalWeight[i])
              sampleMoisture <- as.numeric(finalResult[i])
              dryProfile <- dryingProfile[i]
              modelAbb <- corrDF$modelAbbr[p]
              modelName <- corrDF$modelNames[p]
              modelCorr <- round(as.numeric(corrDF$corr[p]), 8)
              parserVer <- paste(parserID, parserVersion, sep = "")
              modelError <- get(corrDF$modelAbbr[p])[[5]]
              
              modelResultsNew <- data.frame(todaysdate, logTime, datafile = infile, sampleID, dryProfile, dryingTime, initialTemp, finalTemp, sampleMass, sampleMoisture, modelAbb, modelName,
                                            modelCorr, modelEquation,
                                            modelEstimates, modelStdErrors, modeltValues, modelPr.gt.ts, parserVer, modelError)

              if (is.null(modelResults)) {
                modelResults <- modelResultsNew
              } else {
                modelResults <- rbind(modelResults, modelResultsNew)
              }
              
              HTML(paste("Pearson's correlation: ", corrDF$corr[p]), file = HTMLoutput)
              HTML("<hr>", file = HTMLoutput)                
            }
            
            HTML(paste(corrDF), file = HTMLoutput, digits = 8, row.names = FALSE)
            HTML("<hr>", file = HTMLoutput)
            
            # Compile data in .csv
            
            
            sumDataNext <- cbind(rawData, predictDFW, datafile = infile)
            
            if (is.null(sumData)) {
              
              sumData <- sumDataNext
            } else {
              
              # Append next data set from current file to sumData dataframe
              
              sumData <- rbind(sumData, sumDataNext)
            }
            
            # Append raw data to predictDF for export
            
            predictDF <- rbind(plotData, predictDF)
            
            # Include data path
            
            predictDF <- cbind(predictDF, datafile = infile)
            
            # Export model results for topModels
            
            modelOutfile <- paste(dataFolder, modelcsv, sep = "/")
            write.csv(predictDF, file = modelOutfile, quote = TRUE, row.names = FALSE)
            
            # Output sample info to _info.txt
            
            cat(infile, testID[i], switchOffMode[i], dryingProfile[i], dryingTemp[i], resultUnits[i],
                initialWeight[i], finalWeight[i], elapsedTime[i], finalPanTemp[i], finalResult[i],i, sep = ",", file = infoName, append = TRUE)
            cat("\n", file = infoName, append = TRUE)
            
            # Save summary data
            
            outfile <- paste(dataFolder, csvName, sep = "/")
            modFile <- paste(dataFolder, modName, sep = "/")

            if (!is.null(sumData)) {
              write.csv(sumData, file = outfile, quote = TRUE, row.names = FALSE)
            }
            
            # Update parser log
            
            cat(infile, ", ", file = logfName, append = TRUE)
            cat(outfile, "\n", file = logfName, append = TRUE)
            cat("\n\nMoistureBalance parsed data written to:\n\n", outfile, "\n")
            
            HTML.title(paste("<a href=file:///", gsub(" ", "%20", modelOutfile) , ">", modelOutfile, "</a>", sep = ""), HR = 3, file = HTMLoutput)
            HTML.title(paste("<a href=file:///", gsub(" ", "%20", outfile) , ">", outfile, "</a>", sep = ""), HR = 3, file = HTMLoutput)
            #HTMLEndFile()
            write.csv(modelResults, file = modFile, row.names = FALSE)
            
            
            # *** This is where sumData, predictDFW and predictDF should be updated
            
            
            
          } else {
            
            # There are less than 2 data points in the data file
            
            cat(infile, ", Hey, Requires 2 or more data points to continue.\n", file = logfName, append = TRUE)
          }
          
        }
        } else {
          
          # There appears to be a data mismatch in the file
          
          cat(infile, ", Hey, Data mismatch error.\n", file = logfName, append = TRUE)
        }
        
      } else {
        
        # The file contains no data
        
        cat(infile, ", Hey, The file has an unexpected format (not MB45 data).\n", file = logfName, append = TRUE)
      }
      
    } else {
      
      # The file has size less or equal to zero
      
      cat(infile, ", Hey, The file size is less than or equal to zero.\n", file = logfName, append = TRUE)
    }
  
  # Concatenate model results for all datafiles
  
  if (is.null(allModelParams)) {
     allModelParams <- modelResults
  } else {
    allModelParams <- rbind(allModelParams, modelResults)
  }
  
  if (is.null(allModelFits)) {
    allModelFits <- sumData
  } else {
    allModelFits <- rbind(allModelFits, sumData)
  }
  
  if (is.null(allModelFitsLong)) {
    allModelFitsLong <- predictDF
  } else {
    allModelFitsLong <- rbind(allModelFitsLong, predictDF)
  }
  
  # *** Should predictDF and modelResults also be set to NULL?
  
  modelResults <- NULL #
  sumData <- NULL
  predictDF <- NULL #
  
 }

# Output model results for all data files.

rootName <- strsplit(dataFolder, "/")[[1]][length(strsplit(dataFolder, "/")[[1]])]
allModParams <- paste(dataFolder, "/allModelParameters_", rootName, ".csv", sep = "")
write.csv(allModelParams, file = allModParams, row.names = FALSE)

allModFits <- paste(dataFolder, "/allModelFits_", rootName, ".csv", sep = "")
write.csv(allModelFits, file = allModFits, row.names = FALSE)

allModFitsLong <- paste(dataFolder, "/allModelFitsLong_", rootName, ".csv", sep = "")
write.csv(allModelFitsLong, file = allModFitsLong, row.names = FALSE)


#write.csv(modelResults, file = modFile, row.names = FALSE)

  #End of data processing
  
  # Cleanup environment
  
  # rm(beginData,chooseMultiple,corr,corrDF,csvName,dataFolder,decimalTime,dectime,dirList,dryingProfile,dryingTemp,elapsedTime,endData,endFlag,Ex,Exp_fit,f,finalPanTemp,
  #    finalResult,finalWeight,fitData,fName,g,HendersonPabis_fit,HP,HTMLoutput, i,infile,infoCount,infoName,initialWeight,j,k,Le,Lewis_fit,lFinal,Li,Linear_fit,Lo,Log_fit,Loga,
  #    Logarithmic_fit,m,mass,MBData,MHP,mMax,modelAbbr,modelcsv,modelNames,modelOutfile,ModifiedHendersonPabis_fit,ModifiedPage_fit,ModifiedPageII_fit,moisture,
  #    MP,MPII,MR,ms,
  #    nModels,nodename,nSamples,OK.but,outfile,p,Pa,Page_fit,parserCounter,parserID,parserVersion,plotData,plotTitle,pngName,predictDF,predictDFW,predictNext,q,
  #    rawData,logfName,logTime,radiobuttondone,rb1,rb2,rbValue,resultUnits,rowPredDFW,rowPredNext,SF,SimplifiedFick_fit,sumData,switchOffMode,testID,time,tMax,
  #    todaysdate,top,tt,Tt,
  #    TwoTerm_fit,WangSingh_fit,WS,xxx.txt)
  
  #End of data processing
  
  # Cleanup environment
  
  # To generate a list of the variables generated by a script do the following:
  
  # 1) noquote(ls())
  # 2) Copy/Paste output from Console
  # 3) enclose text in double quotes and assign to variable a
  # 4) Execute the following commands in order:
  #       b <- gsub("\\[.+?\\]", "", a)  # Remove square bracket line numbering
  #       c <- gsub("\\s+", " ", b)      # Remove extra spaces between variable names
  #       d <- gsub(" ", ",", c)         # Replace spaces with comma
  # 5) rm(...list of variables in d...)