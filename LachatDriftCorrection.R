# Get list of object in the current Global Environment and make copies of any object that matches an object used in the LachatDriftCorrection.R

# To execute type: source("LachatDriftCorrection.R", print.eval = TRUE)



# Load libraries

library(dplyr)

# Set working directory to the location of the script.

mywd <- "C:/Martin/dataprocessing/R/Rrepo/data-parser"

# Store current working directory and restore setting when script is completed.

currentdir <- getwd()

if(!identical(mywd, currentdir)){
  # Set working directory
  setwd(mywd)
  cat("Working directory set to:", getwd(), "\n\n")
}


# Constants

parserVersion = "LachatDriftCorrection.R"
myFilters = Filters # Substitute .ps in original Filters array with .csv
myFilters[3,1] = "Excel csv files (*.csv)"
myFilters[3,2] = "*.csv"
rownames(myFilters)[3] = "csv"

# Set default data directory

defaultdir <- "L:/LACHAT/Spokas/*.csv"

# Choose data file from default directory or navigate to desired directory in the choose.files() dialog.

infile = choose.files(defaultdir, filters = myFilters[c("csv","All"),], caption = "Choose Lachat datafile")

lachatdata = read.csv(infile, header = TRUE, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)

cat("Successfully read data from ", infile, "\n\n")

# Make copy of data frame for processing purposes.

df <- lachatdata

# Set output filenames by adding the suffix _tidy, _messy or _LDC to the input filename. 

tidyOutfile <- paste(strsplit(infile, ".csv"), "_tidy.csv", sep = "")
messyOutfile <- paste(strsplit(infile, ".csv"), "_messy.csv", sep = "")
LDCOutfile <- paste(strsplit(infile, ".csv"), "_LDC.csv", sep = "")


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

dfLDC <- df

# Store copies of tidy, messy and LDC data.

write.csv(tidy, file = tidyOutfile, quote = FALSE, row.names = FALSE)
write.csv(df, file = messyOutfile, quote = FALSE, row.names = FALSE)
write.csv(dfLDC, file = LDCOutfile, quote = FALSE, row.names = FALSE)

cat("Output data to \n\n", tidyOutfile, "\n", messyOutfile, "\n", LDCOutfile, "\n")


# Clean up variables and detach libraries.

rm(analytes, analyteData, currentdir, defaultdir, df, dfLDC, i, infile, lachatdata, LDCOutfile, messyOutfile, myFilters, mywd, parserVersion, sampleInfo, tidy, tidyOutfile)
detach("package:dplyr", unload=TRUE)