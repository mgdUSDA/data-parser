cat("\n Sunflower QNT PAH data parser.  ", "Last edited: February 16, 2015.\n\n")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Parser for Sparky5973 PCTRPT data generated in HPChemStation
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# To execute type: source("SunPAHQNTParser.R", print.eval = TRUE)

# Clear all variables
rm(list = ls())
require(tcltk)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Setup input and output directories

dir1 = "R:/Sparky4/HPCHEM/1/DATA/"
dirOut = "C:/Martin/Sunflower"
dirWork = "C:/Martin/data processing/R"
dirParser = "C:/Martin/data processing"

timestamp = Sys.time()
timestamp = sub(":", "", timestamp)
timestamp = sub(":", "", timestamp)
dirOut = paste(dirOut, "/", timestamp, "/", sep="")
if(is.na(file.info(dirOut)$isdir) == FALSE){
  #Output folder exists.  Do nothing
}else{
  #Output folder does not exist.  Create new folder.
  dir.create(dirOut, recursive = TRUE)
}

# Constants

systemName = "SPKY4"
parserError = "None"
SumTable = NULL
QSumTable = NULL
datafileName = ""

CompoundNames = c("pTerp", "AcP", "Flu", "Phe", "Ant", "Pyr", "BaA", "Chr", "BbFL", "BkFL", "BaP", "Ind", "DBA", "BP")
dataHeader = c("Information from Data File", "File", "Operator", "Acquired", "Sample Name", "Misc Info", "Vial Number", "CurrentMeth")
qntHeader = c("CompNum", "CompName ", "RT", "QIon", "Response", "Concentration", "Units", "PeakType", "Qvalue")
nCompounds = 14

PctTable = matrix(nrow = nCompounds, ncol = 6)
QntTable = matrix(nrow = 1, ncol = nCompounds+1)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Select MSD Data folder to be parsed.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


cat("Imported MSD SIM Quantification Report data from:\n\n")

dataFolder = choose.dir(default = dir1, caption = "Select folder")
dirList = dir(dataFolder)
nFiles = length(dirList)
SampleSumTable = data.frame(nrow = nFiles, ncol = 8 + nCompounds)
QntTable = data.frame(nrow = nCompounds, ncol = 17)

for (j in 1:nFiles){
	if (length(grep(".D", dirList[j])) == 1){
		infile = paste(dataFolder, "\\", dirList[j], "\\PAHQNT.TXT", sep = "")
		if (file.exists(infile)){
			print(infile)
			datafileName = strsplit(infile, "\\\\")[[1]][7]
			PctRdlns=readLines(infile)
			nLines = length(PctRdlns)

			QntRdlns = PctRdlns[gsub(" ", "", substring(PctRdlns,8,36), fixed = TRUE) %in% CompoundNames]
      
			for (i in 1:nCompounds){
			  QntTable[i, 1] = strsplit(PctRdlns[3],": ", fixed = TRUE)[[1]][2] #Information from Data File:
			  QntTable[i, 2] = strsplit(PctRdlns[4],": ", fixed = TRUE)[[1]][2] #File:
			  QntTable[i, 3] = strsplit(PctRdlns[6],": ", fixed = TRUE)[[1]][2] #Operator:
			  QntTable[i, 4] = strsplit(PctRdlns[5],": ", fixed = TRUE)[[1]][2] #Acquired:
			  QntTable[i, 5] = strsplit(PctRdlns[7],": ", fixed = TRUE)[[1]][2] #Sample Name:
			  QntTable[i, 6] = strsplit(PctRdlns[8],": ", fixed = TRUE)[[1]][2] #Misc Info:
			  QntTable[i, 7] = substring(strsplit(PctRdlns[9],": ", fixed = TRUE)[[1]][2], 1, 4) #Vial Number:
			  QntTable[i, 8] = strsplit(PctRdlns[12],": ", fixed = TRUE)[[1]][2] #CurrentMeth:
			  
			  QntTable[i, 9] = substring(QntRdlns[i],0,6)	#Compound#
			  QntTable[i, 10] = substring(QntRdlns[i],8,36)	#Compound Name 
			  QntTable[i, 11] = substring(QntRdlns[i],37,42)	#Retention Time
			  QntTable[i, 12] = substring(QntRdlns[i],43,46)	#TIC_QIon
			  QntTable[i, 13] = substring(QntRdlns[i],47,55)	#Response
			  QntTable[i, 14] = substring(QntRdlns[i],56,64)	#Concentration   
			  QntTable[i, 15] = substring(QntRdlns[i],65,68)	#Units
			  QntTable[i, 16] = substring(QntRdlns[i],69,71)	#Peak Type   
			  QntTable[i, 17] = substring(QntRdlns[i],72,77)	#Qvalue
			}
      if (!exists("QntSumTable")) {
        QntSumTable = QntTable
      } else {
        QntSumTable = rbind(QntSumTable, QntTable)
      }
      
			SampleSumTable[j, 1] = strsplit(PctRdlns[3],": ", fixed = TRUE)[[1]][2] #jnformation from Data File:
			SampleSumTable[j, 2] = strsplit(PctRdlns[4],": ", fixed = TRUE)[[1]][2] #File:
			SampleSumTable[j, 3] = strsplit(PctRdlns[6],": ", fixed = TRUE)[[1]][2] #Operator:
			SampleSumTable[j, 4] = strsplit(PctRdlns[5],": ", fixed = TRUE)[[1]][2] #Acquired:
			SampleSumTable[j, 5] = strsplit(PctRdlns[7],": ", fixed = TRUE)[[1]][2] #Sample Name:
			SampleSumTable[j, 6] = strsplit(PctRdlns[8],": ", fixed = TRUE)[[1]][2] #Misc Info:
			SampleSumTable[j, 7] = substring(strsplit(PctRdlns[9],": ", fixed = TRUE)[[1]][2], 1, 4) #Vial Number:
			SampleSumTable[j, 8] = strsplit(PctRdlns[12],": ", fixed = TRUE)[[1]][2] #CurrentMeth:
			
			for (i in 1:nCompounds){
			  SampleSumTable[j, i+8] = substring(QntRdlns[i],56,64)
			}
	  }
  }
}

SampleHeader = c(dataHeader, CompoundNames)
colnames(SampleSumTable) = c(dataHeader, CompoundNames)
colnames(QntSumTable) = c(dataHeader, qntHeader)

# outfile = choose.files(dirOut, filters = Filters[c("txt","All"),], caption = "Choose SPARKY 4 output file")

outfile1 = paste(dirOut, "QntSumRes.csv", sep = "")
write.csv(QntSumTable, file = outfile1, quote = TRUE, row.names = FALSE)

outfile2 = paste(dirOut, "SampleSumRes.csv", sep = "")
write.csv(SampleSumTable, file = outfile2, quote = TRUE, row.names = FALSE)

cat("\nSunflower parsed data written to: \n\n", outfile1, "\n\n", outfile2, "\n\n")
tkmessageBox(title = "Sunflower PAH output", message=paste(outfile1, "and", outfile2))
setwd(dirWork)
