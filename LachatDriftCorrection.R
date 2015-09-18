# Get list of object in the current Global Environment and make copies of any object that matches an object used in the LachatDriftCorrection.R

# To execute type: source("LachatDriftCorrection.R", print.eval = TRUE)

# Load libraries

library(dplyr)
library(ggplot2)
library(R2HTML)
library(tcltk)


# Constants
rm(list = ls())

driftThreshold = 0.1 # 10% drift
parserVersion = "LachatDriftCorrection.R"
myFilters = Filters # Substitute .ps in original Filters array with .csv
myFilters[3,1] = "Excel csv files (*.csv)"
myFilters[3,2] = "*.csv"
rownames(myFilters)[3] = "csv"
radiobuttondone <- tclVar(0)
driftCorrect <- 2

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function chooseUnits allows user to choose between concentration and
# mass units for GHG amounts.  ppm are preferred for gas samples; ug for
# insitu incubations.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

doDriftCorrection = function(){
  rbVal = as.numeric(tclvalue(rbValue))
  tkdestroy(tt)
  if (rbVal==1){
    # No longer does anything.  Used to display these cutsie messages...      
    # tkmessageBox(title = "Calibration mode:", message="Auto calibration?  Good luck!")
  }
  if (rbVal==2){
    # tkmessageBox(title = "Calibration mode:", message="Manual calibration? Sorry about that.")
  }
  tclvalue(radiobuttondone) <- 1
  driftCorrect <<- rbVal
}# End of function chooseUnits


# Set working directory to the location of the script.

mywd <- "C:/Martin/dataprocessing/R/Rrepo/data-parser"

# Store current working directory and restore setting when script is completed.


currentdir <- getwd()

if(!identical(mywd, currentdir)){
  # Set working directory
  setwd(mywd)
  cat("Working directory set to:", getwd(), "\n\n")
}


# Set default data directory

defaultdir <- "L:/LACHAT/Spokas/*.csv"

# Choose data file from default directory or navigate to desired directory in the choose.files() dialog.

infile = choose.files(defaultdir, filters = myFilters[c("csv","All"),], caption = "Choose Lachat datafile")

lachatdata = read.csv(infile, header = TRUE, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)

# Only select rows that contain a valid Detection.Date

lachatdata <- lachatdata[!(lachatdata$Detection.Date ==""),] 

cat("Successfully read data from ", infile, "\n\n")

# Make copy of data frame for processing purposes.

df <- lachatdata

# Set output filenames by adding the suffix _tidy, _messy or _LDC to the input filename. 

basefilepath <- strsplit(infile, ".csv")
filename <- strsplit(infile, "\\\\")[[1]][length(strsplit(infile, "\\\\")[[1]])]

tidyOutfile <- paste(basefilepath, "_tidy.csv", sep = "")
messyOutfile <- paste(basefilepath, "_messy.csv", sep = "")
LDCOutfile <- paste(basefilepath, "_LDC.csv", sep = "")

HTMLoutput <- paste(basefilepath, "_report.html", sep = "")
HTML.title(paste("Lachat Drift Correction Report for ", filename, sep=""), Align = "center", HR=3, file=HTMLoutput)
HTML.title(Sys.time(), Align = "center", HR=4, file=HTMLoutput)

# Column names for raw lachat data:
#
#  [1] "Sample.ID"             "Replicate.Number"      "Cup.Number"            "Detection.Date"        "Detection.Time"       
#  [6] "Auto.Dilution.Factor"  "Analyte.Name"          "Peak.Area"             "Peak.Height"           "Peak.Concentration"   
# [11] "Concentration.Units"   "Channel.Number"        "Analyte.Name.1"        "Peak.Area.1"           "Peak.Height.1"        
# [16] "Peak.Concentration.1"  "Concentration.Units.1" "Channel.Number.1"      "Analyte.Name.2"        "Peak.Area.2"          
# [21] "Peak.Height.2"         "Peak.Concentration.2"  "Concentration.Units.2" "Channel.Number.2"


# Identify unique analytes in data file.

analytes <- df %>%
  select(contains("Analyte")) %>%
  unique() %>%
  as.character()

# Rename columns to include analyte name.

colnames(df)[grep("Peak.Area", colnames(df))] <- c(paste(analytes, "Area", sep = "" ))
colnames(df)[grep("Peak.Height", colnames(df))] <- c(paste(analytes, "Height", sep = "" ))
colnames(df)[grep("Peak.Concentration", colnames(df))] <- c(paste(analytes, "Concentration", sep = "" ))
colnames(df)[grep("Concentration.Units", colnames(df))] <- c(paste(analytes, "Units", sep = "" ))
colnames(df)[grep("Channel.Number", colnames(df))] <- c(paste(analytes, "Channel", sep = "" ))

# Remove "." from column names.

colnames(df) <- gsub("[.]", "", colnames(df))

# Manually reshape data.  May be able to do this directly with tidyr::gather().

# Gather sample info.

sampleInfo <- select(df, SampleID:AutoDilutionFactor)

# Gather analyte info and rbind to new data frame
for (i in 1:length(analytes)){
  analyteData <- cbind(select(df, contains(analytes[i])), analytes[i])
  colnames(analyteData) <- c("Area", "Height", "Concentration", "Units", "Channel", "Analyte")
  analyteData <- cbind(sampleInfo, analyteData)
  if (i == 1) {
    # Create new data frame for tidy data
    tidy <- analyteData
  } else {
    tidy <- rbind(tidy, analyteData)  
  }
}


# Correct analyte concentrations for sensitivity variations.  NO3 and PO4 may exhibit some peak area-height drift due to contamination in reagents or manifold lines.  Areas for the highest standard are adjusted as a function of run position.  Sample areas are adjusted according to a spline or polynomial fit of the standard drift.

dfLDC <- tidy

# Convert SampleID to lower case and remove all blank spaces

dfLDC$id <- tolower(dfLDC$SampleID)
dfLDC$id <- gsub("\\s", "", dfLDC$id) 

# Select 20ppm standards which will be used for drift correction if needed. 


Nitrate <- dfLDC[dfLDC$Analyte == "Nitrate",]
Nitrate$index <- as.numeric(row.names(Nitrate))
Nitrate20ppm <- filter(Nitrate, grepl("S", CupNumber) & grepl("20", id) & !grepl("po4", id))


# There need to be at least 2 standards in order to perform the drift correction.

if (nrow(Nitrate20ppm) > 1) {
  # Calculate standard error for the standards

  stdErr20ppm <- mean(Nitrate20ppm$Concentration, na.rm = TRUE) / Nitrate20ppm$Concentration[1]
  
  # Calculate polynomial order
  
  polyOrder <- min(nrow(Nitrate20ppm) - 1, 3)
  
  # Check if 20ppm standard areas change more than threshold value driftThreshold.
  
  cat("Drift correction threshold is", driftThreshold, "\n\n", "20ppm standard error is", stdErr20ppm, "\n\n") 
  
  
  if (driftThreshold < stdErr20ppm) {
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    # Choose whether to perform drift correction or not
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    
    driftCorrect = NA
    tt <- tktoplevel()
    tktitle(tt) <- "Perform drift correction?"
    rb1 <- tkradiobutton(tt)
    rb2 <- tkradiobutton(tt)
    rbValue <- tclVar(1)
    tkconfigure(rb1,variable=rbValue,value=1)
    tkconfigure(rb2,variable=rbValue,value=2)
    tkgrid(tklabel(tt,text="Perform drift correction:"))
    tkgrid(tklabel(tt,text="Yes"),rb1)
    tkgrid(tklabel(tt,text="No"),rb2)
    OK.but <- tkbutton(tt, text="OK", command = doDriftCorrection)
    tkgrid(OK.but)
    tkfocus(tt)
    tkwait.variable(radiobuttondone)
    
    if (driftCorrect == 1) {
      polyFit <- lm(Concentration ~ poly(index, polyOrder), data = Nitrate20ppm)
      windows()
      ggplot(Nitrate20ppm, aes(index, Concentration)) +
        geom_smooth(method = "lm", formula = y ~ poly(x, polyOrder), se = FALSE) +
        geom_point() +
        coord_cartesian(ylim = c(min(Nitrate$Concentration), max(Nitrate$Concentration) * 1.1))
      graphpath <- paste(basefilepath, ".png", sep = "" )
      ggsave(graphpath)
      pngName <- strsplit(graphpath, "\\\\")[[1]][length(strsplit(graphpath, "\\\\")[[1]])]
      dev.off()

      HTML("<hr>",file=HTMLoutput)
      HTMLInsertGraph(pngName, file=HTMLoutput, GraphBorder = 3, Align = "center")
      HTML(summary.lm(polyFit), file=HTMLoutput)
      HTML("Fitted standard values:", file = HTMLoutput)
      
      predFit <- predict(polyFit, data.frame(index = Nitrate$index)) # predict() requires new data be passed as a data frame with matching column name
      correctionFactor <- Nitrate20ppm$Concentration[1] / predFit
      correctionFactor[1] <- 1.0
      Nitrate$LDC <- Nitrate$Concentration * correctionFactor
      
      # Set negative concentrations to zero
      
      Nitrate[Nitrate$LDC < 0, "LDC"] <- 0.0
      dfLDC[row.names(Nitrate), "Concentration"] <- Nitrate$LDC
      
      # Output drift corrected data
      
      write.csv(dfLDC, file = LDCOutfile, quote = FALSE, row.names = FALSE)
    } else {
      cat("Standard Error of 20ppm stds is ", stdErr20ppm, "but, no drift correction performed.\n\n")
    }
  } else {
    cat("Standard Error of 20ppm stds is ", stdErr20ppm, ". No drift correction necessary.\n\n")
  }
  
} else {
  cat("There were less than 2 standards required for drift correction.  No drift correction performed.\n\n")
}
# Store copies of tidy, messy and LDC data.

write.csv(tidy, file = tidyOutfile, quote = FALSE, row.names = FALSE)
write.csv(df, file = messyOutfile, quote = FALSE, row.names = FALSE)


cat("Output data to \n\n", tidyOutfile, "\n", messyOutfile, "\n")
if (driftCorrect == 1) {
  cat(LDCOutfile, "\n")
#  rm(correctionFactor, polyFit, predFit, analytes, analyteData, currentdir, defaultdir, df, dfLDC, driftCorrect, doDriftCorrection, driftThreshold, i, infile, lachatdata, LDCOutfile, messyOutfile, myFilters, mywd, Nitrate, Nitrate20ppm, OK.but, parserVersion, polyOrder, radiobuttondone, rb1, rb2, rbValue, sampleInfo, tidy, tidyOutfile, tt)
} else {
#  rm(analytes, analyteData, currentdir, defaultdir, df, dfLDC, driftCorrect, doDriftCorrection, driftThreshold, i, infile, lachatdata, LDCOutfile, messyOutfile, myFilters, mywd, Nitrate, Nitrate20ppm, parserVersion, radiobuttondone, sampleInfo, tidy, tidyOutfile)
} 
  
# detach("package:dplyr", unload = TRUE)
# detach("package:ggplot2", unload = TRUE)
# detach("package:R2HTML", unload = TRUE)
# detach("package:tcltk", unload = TRUE)