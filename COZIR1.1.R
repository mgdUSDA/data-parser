cat("\n COZIR1.1 Last edited: September 17, 2015.\n\n")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Parses Sartorius balance data and creates a plot
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Critical problems:
#
# sampleTime interval checking.  Currently works for sample interval 10 seconds
#  
# 

# Execute ---------------------------------------------------------------------
# To execute type: source("COZIR1.1.R", print.eval = TRUE)

# Clear variables and console

rm(list = ls())
cat("\014")

#library(dplyr)
library(ggplot2)
library(R2HTML)
library(tcltk)


parserVersion = "COZIR September 17, 2015."

# Function definition ----


getSampleTime = function(){
  # get date/time for each sample reading rounded to the nearest second
  # based on interval.  
  
  st <- as.POSIXct(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  st.Sec <- as.numeric(format(st, "%S"))
  while(st.Sec %% interval != 0){
    st <- as.POSIXct(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    st.Sec <- as.numeric(format(st, "%S"))
    #  print(paste("Wait for it...[", st, "]"))
  }
  return(st)
}


graphData <- function(cozirdata, path){
  graphTitle <- strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])]
  windows()
  ggplot(cozirdata, aes(sampletime, co2)) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1, se = FALSE) +
    geom_point() +
    ggtitle(paste("CO2 data for ", graphTitle, sep = "")) +
    ylab("CO2 concentration (ppm)") +
    xlab("time (min)")
  ggsave(path)
  dev.off()
}


getFlux <- function(cozirdata, path) {
  graphTitle <- strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])]
  lmFit <- lm(co2 ~ sampletime, data = cozirdata)
  xmax <- round(max(cozirdata$sampletime, na.rm = TRUE))
  ymax <- round(max(cozirdata$co2, na.rm = TRUE))
  xpos <- round(0.3 * xmax)
  ypos <- round(0.9 * ymax)
  windows()
  ggplot(cozirdata, aes(sampletime, co2)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#990000", size = 1, se = FALSE) +
    geom_point() +
    ggtitle(paste("CO2 flux for", graphTitle)) +
    ylab("CO2 concentration (ppm)") +
    xlab("time (min)") +
    annotate("text", x = xpos, y = ypos, label = paste("flux =",
              round(lmFit[[1]][2], digits = 1), "ppm CO2/min"))
  ggsave(path)
  dev.off()
  return(lmFit[[1]][2])
}

# Resets variable values when dialog button is pressed

reset <- function(){
  tclvalue(site) <- "COZIR"
  tclvalue(interval) <- "10"
  tclvalue(duration) <- "600"
  tclvalue(nPorts) <- "1"
}


# Submits variable values when dialog button is pressed

submit <- function() {
  m0 <- tclvalue(site)
  m1 <- as.numeric(tclvalue(interval))
  m2 <- as.numeric(tclvalue(duration))
  m3 <- as.numeric(tclvalue(nPorts))
  
  #  tkmessageBox(message="Done!")
  tkdestroy(tt)
  tclvalue(didone) = 1
  site <<- m0
  interval <<- m1
  duration <<- m2
  nPorts <<- m3
}


# End of function definitions

# Main ----
site = "COZIR"
duration = 600
interval = 10


didone <- tclVar(0)
nPorts <- 1 # Number of sensors/ports

# Setup data output directory and file name: dirOut() and fname

reportTime <- Sys.time()
reportDate <- substr(reportTime, 1, 10)
reportTime <- sub(":", "", reportTime)
reportTime <- sub(":", "", reportTime) # Once for each ":" in the time.
reportTime <- strsplit(reportTime, " ")[[1]][2]

dirOut = paste("C:/COZIR", reportDate, sep="/")
if(is.na(file.info(dirOut)$isdir) == FALSE){
  #Output folder exists.  Do nothing
}else{
  #Output folder does not exist.  Create new folder.
  dir.create(dirOut, recursive = TRUE)
}


newData <- data.frame(matrix(nrow = 1, ncol = 6))
colnames(newData) <- c("datetime", "sampletime", "temperature", "humidity", "co2", "o2")


system("MODE COM8:9600, N, 8, 1")
port1 = file("COM8", open = 'r+')
isOpen(port1)


# Get sampling interval (s) and number of sensors.----
site <- tclVar("COZIR")
interval <- tclVar("10")
duration <- tclVar("600")
nPorts <- tclVar("1")

tt <- tktoplevel()
tkwm.title(tt,"Soil mass input")
m0.entry <- tkentry(tt, textvariable = site)
m1.entry <- tkentry(tt, textvariable = interval)
m2.entry <- tkentry(tt, textvariable = duration)
m3.entry <- tkentry(tt, textvariable = nPorts)
reset.but <- tkbutton(tt, text = "Reset", command = reset)
submit.but <- tkbutton(tt, text = "Submit", command = submit)
tkgrid(tklabel(tt, text = "Reading setup:"),columnspan=2)
tkgrid(tklabel(tt, text = "Site name"), m0.entry)
tkgrid(tklabel(tt, text = "Sample interval (s)"), m1.entry)
tkgrid(tklabel(tt, text = "Sample duration (s)"), m2.entry)
tkgrid(tklabel(tt, text = "Number of sensors (max. 2)"), m3.entry)
tkgrid(submit.but, reset.but)
tkwait.variable(didone)


# Set data output file paths based on default dirOutC:/COZIR/COZIR_<reportTime>.csv

outfile <- paste(dirOut, "/", site, "_", reportTime, ".csv", sep = "")
HTMLoutput <- paste(strsplit(outfile, "\\.")[[1]][1], ".html", sep = "")
graphPNG <- paste(strsplit(outfile, "\\.")[[1]][1], ".png", sep = "")
fluxPNG <- paste(strsplit(outfile, "\\.")[[1]][1], "_flux.png", sep = "")


# Calculate Number of samples

nSamples = round(duration / interval)

for(i in 1:nSamples){
  sampleTime = getSampleTime()
  dt <- as.character(format(sampleTime, "%Y-%m-%d %H:%M:%S"))
  if(i == 1) {
    dt0 <- dt
  }
  cat("[", i, "] :", dt, "\n")
  newData$datetime <- dt
  newData$sampletime <- as.numeric(difftime(dt, dt0, units = "mins"))
  tryCatch(isOpen(port1),
           warning = function(w) {
             print(paste(summary(port1)$description, "is not open. Reset port."))
             port1 = file("COM8", open = 'r+')
           },
           error = function(e) {
             print(paste(summary(port1)$description, "is not open. Reset port."))
             port1 = file("COM8", open = 'r+')
           }
  )
  
  write("T\r\n", port1) # Return most recent temperature measure
#  Sys.sleep(0.3)
  t <- scan(port1, n = 3, what = "character", quiet = TRUE)
  newData$temperature <- ((as.numeric(t[2])) %% 1000) / 10.0
  cat("T   :", newData$temperature, "\n")

  write("H\r\n", port1) # Return most recent humidity measure
#  Sys.sleep(0.3)
  h <- scan(port1, n = 3, what = "character", quiet = TRUE)
  newData$humidity <- as.numeric(h[2]) / 10.0
  cat("H   :", newData$humidity, "\n")
  
  write("Z\r\n", port1) # Return most recent CO2 measure
#  Sys.sleep(0.3)
  c <- scan(port1, n = 3, what = "character", quiet = TRUE)
  newData$co2 <- as.numeric(c[2])
  cat("CO2 :", newData$co2, "\n")
  
# repeat port check for O2 sensor.  
  if (i == 1){
    compData <- newData[1,]
  }else{
    compData <- rbind(compData, newData)
  }

Sys.sleep(1)

}
close(port1)
write.csv(compData, file = outfile, quote = FALSE, row.names = FALSE)
graphData(compData, graphPNG)
newFlux <- getFlux(compData, fluxPNG)

HTML.title(paste("COZIR CO2 sensor reading report for ", site, "_", reportTime, sep=""),
           Align = "center", HR=3, file=HTMLoutput)
HTML.title(reportDate, Align = "center", HR=4, file=HTMLoutput)
HTML("<hr>",file=HTMLoutput)
HTML.title(paste(site, "flux = ", round(newFlux, digits = 1), "ppm CO2/min."),
           Align = "center", HR=4, file=HTMLoutput)

pngName <- strsplit(graphPNG, "/")[[1]][length(strsplit(graphPNG, "/")[[1]])]
HTMLInsertGraph(pngName, file=HTMLoutput, GraphBorder = 3, Align = "center")
pngName <- strsplit(fluxPNG, "/")[[1]][length(strsplit(fluxPNG, "/")[[1]])]
HTMLInsertGraph(pngName, file=HTMLoutput, GraphBorder = 3, Align = "center")

browseURL(HTMLoutput)

# Cleanup Global Environment variables and libraries

# rm("c", "compData", "didone", "dirOut", "dt", "dt0", "duration", "getFlux", "getSampleTime",
#    "graphData",  "h", "i", "interval", "m0.entry",  "m1.entry",  "m2.entry", "m3.entry",
#    "newData", "newFlux",  "nPorts",  "nSamples", "outfile", "parserVersion", "path",
#    "pngName", "port1", "reportDate", "reportTime",  "reset", "reset.but", "sampleTime",
#    "site", "submit", "submit.but", "t", "tt")
# 
# detach("package:dplyr", unload = TRUE)
# detach("package:ggplot2", unload = TRUE)
# detach("package:R2HTML", unload = TRUE)
# detach("package:tcltk", unload = TRUE)
