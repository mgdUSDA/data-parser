cat("\n BlackBear data parser.  ", "Last edited: June 25, 2015.\n\n")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Parser for Sparky4 PCTRPT data generated in HPChemStation
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# To execute type: source("BlackBearParser15.6.R", print.eval = TRUE)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# To do list ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Process all files in a directory

# Report whether data file contains peaks or not.

# Clear all variables
rm(list = ls())
require(tcltk)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Setup input and output directories
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dir1 = "C:/Martin/BlackBear/HPCHEM/1/DATA/"
dir2 = "C:/Martin/BlackBear/HPCHEM/1/DATA/"
dirOut = "C:/Labwork/BlackBear/"
dirWork = "C:/Martin/data processing/R/"
dirParser = "C:/Martin/data processing/"


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Constants
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

parserVersion = "BlackBear15.6.1.R"
systemName = "BlackBear"
parserError = "None"
SumTable = NULL
SumData = NULL
SumResults = NULL
nChar_RefNum = 6
nChar_CAS = 11
nChar_Qual = 3

round2 = function(x, n) {
  z = NULL
  z_frac = NULL
  z = abs(x)*10^n
  z
  z_frac = z - trunc(z)
  z_frac
  z_frac = round(z_frac, digits = 2)
  for (i in 1:length(z)){
    if (z_frac[i] >= 0.5) { # Couldn't get if condition to apply to each element of z_frac
      z[i] = trunc(z)[i] + 1
    }else{
      z[i] = trunc(z)[i]
    }
    z[i] = z[i]/10^n
    z[i]*sign(x[i])
  }
  return(z)
}

trimWhiteSpace = function(x) {
  x = sub("^\\s+", "", x) # trim leading white spaces
  x = sub("\\s+$", "", x) # trim trailing white spaces
  return(x)  
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The maximum number of library matches reported per peak is set in HPChem = 3.
# nLibResults is the number of results reported by the parser in LibData.  Sample summary
# will only report 1 compound per peak.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

nLibResults = 2


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Select MSD data folder to be processesed.
# The folder must contain a LIBRPT file and may also contain
# a PCTRPT file
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cat("Imported MSD Library Report data from:\n\n")

dataFolder = choose.dir(default = dir1, caption = "Select folder")
seqName = strsplit(dataFolder, "\\\\")[[1]][length(strsplit(dataFolder, "\\\\")[[1]])]

reportTime = Sys.time()
reportTime = sub(":", "", reportTime)
reportTime = sub(":", "", reportTime) # Once for each ":" in the time.
dirOut = paste(dirOut, seqName, "/", reportTime, "/", sep="")
if(is.na(file.info(dirOut)$isdir) == FALSE){
  #Output folder exists.  Do nothing
}else{
  #Output folder does not exist.  Create new folder.
  dir.create(dirOut, recursive = TRUE)
}
repfName = paste(dirOut, "VOCReport_", reportTime, ".txt", sep = "")
cat("Parser version:", parserVersion, "\n", file = repfName, append = TRUE)
cat("Sequence name:", dataFolder, "\n", file = repfName, append = TRUE)
cat("Date:", date(), "\n", file = repfName, append = TRUE)

dirList = dir(dataFolder)
nFiles = length(dirList)

xxx.D = grep(".D", dirList)
for (j in xxx.D){
  cat("Processing data for file", j, ":", dirList[j], "\n")
  
  infileLib = paste(dataFolder, "\\", dirList[j], "\\LIBRPT.TXT", sep = "")
  infilePct = paste(dataFolder, "\\", dirList[j], "\\PCTRPT.TXT", sep = "")
  
  # Check if both files exist.  If only LIBRPT exists, process data and fill PCTRPT gaps with NA.
  
  if (file.exists(infileLib)){ # If LIBRPT exists
    
    cat(dirList[j], ": LIBRPT.TXT processed.\n", file = repfName, append = TRUE)
    cat(dirList[j], ": LIBRPT.TXT processed.\n")
    print(infileLib)
    Rdlns=readLines(infileLib)
    nLines = length(Rdlns)
    
    # Find first row of data. In Library Search report, a hard rule is printed directly 
    # above the first row of data.  Use grep to detect that hard rule.
    
    if (length(grep("________", Rdlns)) == 0) {
      cat(dirList[j], ": No peaks detected. Skip PCTRPT.\n", file = repfName, append = TRUE)
      cat(dirList[j], ": No peaks detected. Skip PCTRPT.\n")
      NoPeaks = TRUE
      
      # Insert a blank row in the summary table
      
      # Assume that if LIBRPT is empty so is PCTRPT
    } else {
      startRowLib = min(grep("________", Rdlns))  
      
      # Copy sample info: LIBRPT header plus the last line of the file.
      
      LibInfo = Rdlns[c(1:(startRowLib-2), nLines)]
      LibRdlns = Rdlns[(startRowLib+1):(nLines-1)]
      nDataLines = length(LibRdlns)
      nrowsLib = nDataLines
      LibTable = matrix(nrow = nDataLines, ncol = 7)
      colnames(LibTable) = c("PkNum", "RT", "AreaPct", "LibraryID", "RefNum", "CASNum", "Qual")
      for (i in 1:nDataLines){
        if (nchar(LibRdlns[i]) > 10){
          if (length(grep("DATABASE", LibRdlns[i])) == 1){
            PkNum = substring(LibRdlns[i], 1, 3)   	  # Pk number
            RT = substring(LibRdlns[i], 5, 10)		    # RT
            AreaPct = substring(LibRdlns[i], 12, 16)	# Area%
          }else{
            LibTable[i, 1] = PkNum	                         # Pk number repeated
            LibTable[i, 2] = RT	                             # RT repeated
            LibTable[i, 3] = AreaPct	                       # Area% repeated
            if (nchar(LibRdlns[i]) == 73){ # Library/ID does NOT extends past column 51
              LibTable[i, 4] = substring(LibRdlns[i], 18, 51)  # Library/ID
              LibTable[i, 5] = substring(LibRdlns[i], 53, 58)  # Reference number
              LibTable[i, 6] = substring(LibRdlns[i], 60, 70)  # CAS number
              LibTable[i, 7] = substring(LibRdlns[i], 72, 73)  # Quality
            }else{
              indexShift = nchar(LibRdlns[i]) - 73
              LibTable[i, 4] =  substring(LibRdlns[i], 18, 51 + indexShift)
              LibTable[i, 5] = substring(LibRdlns[i], 53 + indexShift, 58 + indexShift)  # Reference number
              LibTable[i, 6] = substring(LibRdlns[i], 60 + indexShift, 70 + indexShift)  # CAS number
              LibTable[i, 7] = substring(LibRdlns[i], 72 + indexShift, 73 + indexShift)  # Quality
            }
          }
        }
      }
      
      # Ignore blank lines and copy as data frame 
      
      LibData = data.frame(LibTable[!is.na(LibTable[,1]),], stringsAsFactors = FALSE)
      
      # Convert factor to numeric
      
      LibData$PkNum = as.numeric(as.character(LibData$PkNum))
      LibData$RT = as.numeric(as.character(LibData$RT))
      LibData$AreaPct = as.numeric(as.character(LibData$AreaPct))
      LibData$RefNum = as.numeric(as.character(LibData$RefNum))
      LibData$Qual = as.numeric(as.character(LibData$Qual))
      
      # Subset only top compound per peak by counting how many compounds are
      # reported per peak, N and then skipping every N-th row.
      
      uPkNum = unique(LibData$PkNum)
      topResults = sapply(uPkNum, function(x) min(which(LibData$PkNum == x)))
      topLibData = LibData[topResults,]
      
      # Look for PCTRPT.TXT only if LIBRPT.TXT exists AND NoPeaks is FALSE
      
      if (file.exists(infilePct)){ # If PCTRPT exists
        
        cat(dirList[j], ": PCTRPT.TXT processed.\n", file = repfName, append = TRUE)
        cat(dirList[j], ": PCTRPT.TXT processed.\n")
        print(infilePct)
        Rdlns = readLines(infilePct)
        nLines = length(Rdlns)
        minRowPct = min(grep("_______", Rdlns))
        maxRowPct = max(grep("_______", Rdlns))
        startRowPct =  (grep("_______", Rdlns))[2]
        PctInfo = Rdlns[c((minRowPct+1):(startRowPct-1), nLines)]
        PctRdlns = Rdlns[(startRowPct+3):(maxRowPct-1)]
        
        # Extract sample info and include in SumData
        
        sampleName = strsplit(PctInfo[grep("Sample Name", PctInfo)], ":")[[1]][2]
        sampleName = trimWhiteSpace(sampleName)
        
        fileName = strsplit(PctInfo[grep("File  ", PctInfo)], ":")[[1]][3]
        fileName = trimWhiteSpace(fileName)
        method = strsplit(PctInfo[grep("CurrentMeth", PctInfo)], ":")[[1]][3]
        method = trimWhiteSpace(method)
        vial = strsplit(PctInfo[grep("Vial Number", PctInfo)], ":")[[1]][2]
        vial = as.numeric(as.character(trimWhiteSpace(vial)))
        dateString = PctInfo[length(PctInfo)]
        dateString = trimWhiteSpace(dateString)
        dateTime = strptime(dateString, "%a %b %d %H:%M:%S %Y")
        
        # Now remove hard rules
        
        PctRdlns = PctRdlns[!(seq(1, length(PctRdlns)) %in% grep("________", PctRdlns))]
        PctTable = matrix(nrow = length(PctRdlns), ncol = 5)
        colnames(PctTable) = c("RT", "SignalDescr", "Area", "PkPct", "LPkPct")
        for (i in 1:length(PctRdlns)) {
          PctTable[i, 1] = substring(PctRdlns[i],1,16)   # Ret Time
          PctTable[i, 2] = substring(PctRdlns[i],17,32)   # Signal Descr
          PctTable[i, 3] = substring(PctRdlns[i],33,45)  # Area
          PctTable[i, 4] = substring(PctRdlns[i],46,56)  # %Pk
          PctTable[i, 5] = substring(PctRdlns[i],57,68)	 # %LPk
        }
        
        # Copy as data frame
        
        PctData = data.frame(PctTable, stringsAsFactors = FALSE)
        
        #   # Convert factor RT to numeric and round so that retention times match LibData
        #   # Tried to use a function to round RT but had trouble applying z_frac >= 0.5 
        #   # condition to each element of z_frac
        #   
        #   PctRT = round(as.numeric(as.character(PctData$RT)), digits = 1)
        #   n = 2
        #   z = abs(PctRT)*10^n
        #   
        #   z_frac = z - trunc(z)
        #   for (i in 1:length(z)){
        #     if (z_frac[i] >= 0.5) {
        #       z[i] = trunc(z[i]) + 1
        #     }else{
        #       z[i] = trunc(z[i])
        #     }
        #     z[i] = z[i]/10^n
        #     z[i]*sign(PctRT[i])
        #   }
        #   
        #   PctData$RT = z
        
        # Convert remaining factor to numeric
        
        PctData$RT = as.numeric(as.character(PctData$RT))
        PctData$Area = as.numeric(as.character(PctData$Area))
        PctData$PkPct = as.numeric(as.character(PctData$PkPct))
        PctData$LPkPct = as.numeric(as.character(PctData$LPkPct))
        
        # Merge PctData and topLibData
        
        SumData = cbind(topLibData, PctData)
        
        # Append columns with sample info from PCTRPT.TXT
        
        SumData$dateTime = dateTime
        SumData$vial = vial
        SumData$sampleName = sampleName
        SumData$fileName = fileName
        SumData$method = method
        
        if (is.null(SumResults)) {
          SumResults = SumData
          BlankRow = data.frame(matrix("", nrow = 1, ncol = ncol(SumData)))
          colnames(BlankRow) = names(SumData)
          BlankRow$dateTime[1] = NA
        } else {
          SumResults = rbind(SumResults, SumData, BlankRow)
        }
        
      } else {
        
        cat(dirList[j], ": PCTRPT.TXT missing.\n", file = repfName, append = TRUE)
        cat(dirList[j], ": PCTRPT.TXT missing.\n")
        
      }
      
    } # End of LIBRPT and PCTRPT processing loop
    
    
  } else {
    
    cat(dirList[j], ": LIBRPT.TXT missing.\n", file = repfName, append = TRUE)
    cat(dirList[j], ": LIBRPT.TXT missing.\n")
    
  } # End of loop for processing LIBRPT.TXT
  
  
} # End of directory processing loop

cat("Output summary of composite data.")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Output merged data file to .csv ----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

outfile = choose.files(paste(dirOut, "*.csv", sep = ""), filters = Filters[c("txt","All"),], caption = "Choose BlackBear output file")
write.csv(SumResults, file = outfile, quote = TRUE, row.names = FALSE)

cat("\n\nBlackBear parsed data written to:\n\n", file = repfName, append = TRUE)
cat(outfile, "\n", file = repfName, append = TRUE)
cat("\n\nBlackBear parsed data written to:\n\n", outfile, "\n")

setwd(dirWork)
