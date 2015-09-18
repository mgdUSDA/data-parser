cat("\n Sparky Systems data parser.  ", "Last edited: June 24, 2014.\n\n")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Parses SparkyJr summary data generated in HPChemStation
# and ouputs fit sample data based on least squares linear
# regression.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Critical problems:
# TC:  The parser is generating an extra column in the "long" output file.
# 
# 
# 
# 
#
# Suggested improvements:
# ALL:  Output std names when reporting residuals and expected values
# ALL:  Add sample type column and identify standards with "STD"
# ALL:  Calculate average AC CO2 conc and subtract from fit sample values
# ALL:  Look-up matching FIDTCD or SRI file based on folder name for ECDFID.txt file[requires testing from S526]
# ALL:  Progress bar
# ALL:  Run log:  Sequence table name MMDDYYC, system (JR, S3, TC), date-time, user, completed
# 
# TC:  Implement auto calibration option for TotalChrom stover incubation data.
# TC:  Eliminate NA column at the end of SumTable
# 



# To execute type: source("SPKYNetFit14.4.R", print.eval = TRUE)

# Clear all variables
rm(list = ls())
require(tcltk)
library("R2HTML")
cat("\014")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Constants

setClass(Class="Regress",
         representation(
           lm="lm",
           var="character",
           rsq="numeric"
         )
)

setClass(Class="mydirlist",
         representation(
           d1="character",
           d2="character",
           dO="character",
           dP="character",
           dR="character",
           dW="character"
         )
)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Setup input and output directories
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Network Run from Labs526

# dir1 = "//10.10.10.2/HPCHEM/1/DATA/*.txt"
# dir2 = "//10.10.10.2/HPCHEM/2/DATA/*.txt"
# dirOut = "C:/Users/USDA-ARS-Soils526/LabWork/SPARKY3/*.csv"
# dirWork = "C:/Users/USDA-ARS-Soils526/LabWork/Data processing/R"
# dirParser = "C:/Users/USDA-ARS-Soils526/Labwork/Data processing"


# Constants
parserVersion = "SPKYNetFit14.4.R"
myFilters = Filters # Substitute .ps in original Filters array with .csv
myFilters[3,1] = "Excel csv files (*.csv)"
myFilters[3,2] = "*.csv"
rownames(myFilters)[3] = "csv"

nCompounds = 11
colsPerComp = 8 # Number of columns per compound in HPCHEM reports
dir.new = FALSE
parserError = "None"
pValue = 0.05
computer = "T1700"
recalThresh = 0.05 # NRMSE > 5% triggers recalibration
min.cal = 2 # Minimum number of unique stds used in calibration regression
radiobuttondone <- tclVar(0)

fitCompTable =  matrix(nrow = nCompounds, ncol = 3)
fitCompTable[1,] = c("cal.CO2","Area_CO2", "Height_CO2")
fitCompTable[2,] = c("cal.N2O","Area_N2O", "Height_N2O")
fitCompTable[3,] = c("cal.CH4","Area_CH4", "Height_CH4")
fitCompTable[4,] = c("cal.O2","Area_O2", "Height_O2")
fitCompTable[5,] = c("cal.N2","Area_N2", "Height_N2")
fitCompTable[6,] = c("cal.C2H2","Area_C2H2", "Height_C2H2")
fitCompTable[7,] = c("cal.C2H4","Area_C2H4", "Height_C2H4")
fitCompTable[8,] = c("cal.C2H6","Area_C2H6", "Height_C2H6")
fitCompTable[9,] = c("cal.CO2_ECD","Area_CO2_ECD", "Height_CO2_ECD")
fitCompTable[10,] = c("cal.N2O_ECD","Area_N2O_ECD", "Height_N2O_ECD")
fitCompTable[11,] = c("cal.CH4_PID","Area_N2O_PID", "Height_N2O_PID")

colnames(fitCompTable) = c("cal_name", "Area_name", "Height_name")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Initialize and fill standard arrays with values corresponding to stds list based
# on Sample Name.  For now hardcoding might be the easiest way to do this.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



CO2_ECD_std = c(NA) #Usually only useful for greater than 1% CO2
CO2_std = c(NA)
N2O_ECD_std = c(NA)
N2O_PID_std = c(NA) # Only present at very high concentrations
CH4_std = c(NA)
CH4_PID_std = c(NA) # Only present at very high concentrations
C2H2n_std = c(NA) # Acetylene, ethylene, ethane
O2_std = c(NA)
N2_std = c(NA)

gc = NA
calmode = NA
units = NA




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function definitions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function Two_point.lm() Adapted from flux package [flux.odae.R]
# Two_point.lm() calculates the linear model based on 2-point linear fit.
#
# Example: Two_point.lm("CH4_std", concentrations, area, height, max nrmse = 0.005)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Two_point.lm = function(name_y, conc, a, h, mx.nrms){
  # Put compound data into a data frame for easy manipulation
  gas.df = data.frame(conc, a, h)
  
  # Select all combinations of 2 or more unique stds to generate calibation fits in vers
  vers = unlist(lapply(c(min.cal:length(conc)), function(x) combn(c(1:length(conc)), x, simplify=FALSE)), recursive=FALSE)
  sel = unlist(lapply(vers, function(x) length(unique(gas.df[x, 1])) > 1))
  print(paste("2point sel = ", sum(sel=="TRUE")))
  print(paste("2point sel =", sel))
  vers <- vers[sel] 
  
  # Determine the number of stds per version
  vers.n <- sapply(vers, "length")
  
  # Calculate lm for area and height
  vers.lm.a <- lapply(vers, function(x) lm(conc ~ a, data=gas.df[x,]))
  vers.lm.h <- lapply(vers, function(x) lm(conc ~ h, data=gas.df[x,]))
  
  # Extract the statistic r.squared and calculate the nrmse for lm
  rsq.a = unlist(lapply(vers.lm.a, function(x) summary(x)$r.squared))
  rsq.h = unlist(lapply(vers.lm.h, function(x) summary(x)$r.squared))
  nrmse.a <- unlist(lapply(vers.lm.a, function(x) sqrt(sum(residuals(x)^2)/summary(x)$df[2])/diff(range(x$model[1], na.rm=TRUE))))
  nrmse.h <- unlist(lapply(vers.lm.h, function(x) sqrt(sum(residuals(x)^2)/summary(x)$df[2])/diff(range(x$model[1], na.rm=TRUE))))
  
  # Rank according to nrmse
  # In flux.odae.R,  the lm are ranked first by vers.n then by 1-nrmse favoring the lm that uses more data
  # for example: nrmse.rnk.a <- order(vers.n, 1-nrmse.a, decreasing=TRUE)
  # The calculation for nrmse already takes the range of the fit into account and forcing vers.n may lead to
  # sub-optimum selection of best calibration based on r-squared
  
  # rsq.rnk.a <- order(rsq.a, decreasing=TRUE)
  # rsq.rnk.h <- order(rsq.h, decreasing=TRUE)
  nrmse.rnk.a <- order(1-nrmse.a, decreasing=TRUE)
  nrmse.rnk.h <- order(1-nrmse.h, decreasing=TRUE)
  
  print(paste("nrmse:", nrmse.rnk.a, nrmse.rnk.h))
  
  # Not sure why they used such a complicated calculation for selecting the best calibration:
  # m2t.a = nrmse.rnk.a[nrmse.a[nrmse.rnk.a] <= max.nrmse][1]
  # m2t.h = nrmse.rnk.h[nrmse.h[nrmse.rnk.h] <= max.nrmse][1]
  
  # Select the best calibration for area and height 
  if (nrmse.a[nrmse.rnk.a[1]] <= nrmse.h[nrmse.rnk.h[1]]) {
    m2t = nrmse.rnk.a[1]
    bestVar = "area"
    print(paste("bestCal: Area vrs", name_y))
    std.lm = vers.lm.a[m2t]
    x0 = a
    x = a[unlist(vers[m2t])]
  } else {
    m2t = nrmse.rnk.h[1]
    bestVar = "height"
    print(paste("bestCal: Height vrs", name_y))
    std.lm = vers.lm.h[m2t]
    x0 = h
    x = h[unlist(vers[m2t])]   
  }
  y0 = conc
  y = conc[unlist(vers[m2t])]
  std.lm = lm(y ~ x)
  plotModel(y, x, std.lm, name_y, bestVar, y0, x0)
  return(new("Regress", lm = std.lm, var = bestVar, rsq = summary(std.lm)$r.squared))
} # End of function best.lm()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# This function adds the std values to stdval_cal.XX column for later manual recalibration
# Currently under development AKA not working!!!
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

addStdVal = function(){
  which(colnames(SumTable) == "stdval_cal.CH4")
  SumTable[K5, "Sample Name"]
  CH4_std[c(K5, K3)]
  
} # End of function addStdVal()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function best.lm() Adapted from flux package [flux.odae.R]
# best.lm() calculates the best linear model from all unique combinations of
# of 3 or more standards.
#
# Example: best.lm("CH4_std", concentrations, area, height, max nrmse = 0.005)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

best.lm = function(name_y, conc, a, h, mx.nrms){
  # Put compound data into a data frame for easy manipulation
  gas.df = data.frame(conc, a, h)
  
  # Select all combinations of 2 or more unique stds to generate calibation fits in vers
  vers = unlist(lapply(c(min.cal:length(conc)), function(x) combn(c(1:length(conc)), x, simplify=FALSE)), recursive=FALSE)
  sel = unlist(lapply(vers, function(x) length(unique(gas.df[x, 1])) > min.cal))
  print(paste("sel = ", sum(sel=="TRUE")))
  vers <- vers[sel] 
  
  # Determine the number of stds per version
  vers.n <- sapply(vers, "length")
  
  # Calculate lm for area and height
  vers.lm.a <- lapply(vers, function(x) lm(conc ~ a, data=gas.df[x,]))
  vers.lm.h <- lapply(vers, function(x) lm(conc ~ h, data=gas.df[x,]))
  
  # Extract the statistic r.squared and calculate the nrmse for lm
  rsq.a = unlist(lapply(vers.lm.a, function(x) summary(x)$r.squared))
  rsq.h = unlist(lapply(vers.lm.h, function(x) summary(x)$r.squared))
  nrmse.a <- unlist(lapply(vers.lm.a, function(x) sqrt(sum(residuals(x)^2)/summary(x)$df[2])/diff(range(x$model[1], na.rm=TRUE))))
  nrmse.h <- unlist(lapply(vers.lm.h, function(x) sqrt(sum(residuals(x)^2)/summary(x)$df[2])/diff(range(x$model[1], na.rm=TRUE))))
  
  # Rank according to nrmse
  # In flux.odae.R,  the lm are ranked first by vers.n then by 1-nrmse favoring the lm that uses more data
  # for example: nrmse.rnk.a <- order(vers.n, 1-nrmse.a, decreasing=TRUE)
  # The calculation for nrmse already takes the range of the fit into account and forcing vers.n may lead to
  # sub-optimum selection of best calibration based on r-squared
  
  # rsq.rnk.a <- order(rsq.a, decreasing=TRUE)
  # rsq.rnk.h <- order(rsq.h, decreasing=TRUE)
  nrmse.rnk.a <- order(1-nrmse.a, decreasing=TRUE)
  nrmse.rnk.h <- order(1-nrmse.h, decreasing=TRUE)
  
  print(paste("nrmse:", nrmse.rnk.a, nrmse.rnk.h))
  
  # Not sure why they used such a complicated calculation for selecting the best calibration:
  # m2t.a = nrmse.rnk.a[nrmse.a[nrmse.rnk.a] <= max.nrmse][1]
  # m2t.h = nrmse.rnk.h[nrmse.h[nrmse.rnk.h] <= max.nrmse][1]
  
  # Select the best calibration for area and height 
  if (nrmse.a[nrmse.rnk.a[1]] <= nrmse.h[nrmse.rnk.h[1]]) {
    m2t = nrmse.rnk.a[1]
    bestVar = "area"
    print(paste("bestCal: Area vrs", name_y))
    std.lm = vers.lm.a[m2t]
    x0 = a
    x = a[unlist(vers[m2t])]
  } else {
    m2t = nrmse.rnk.h[1]
    bestVar = "height"
    print(paste("bestCal: Height vrs", name_y))
    std.lm = vers.lm.h[m2t]
    x0 = h
    x = h[unlist(vers[m2t])]   
  }
  y0 = conc
  y = conc[unlist(vers[m2t])]
  std.lm = lm(y ~ x)
  plotModel(y, x, std.lm, name_y, bestVar, y0, x0)
  return(new("Regress", lm = std.lm, var = bestVar, rsq = summary(std.lm)$r.squared))
} # End of function best.lm()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function chooseUnits allows user to choose between concentration and
# mass units for GHG amounts.  ppm are preferred for gas samples; ug for
# insitu incubations.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

chooseCalMode = function(){
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
  calmode <<- rbVal
}# End of function chooseUnits



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function chooseSystems allows user to choose between Sparky Jr. and
# Sparky III
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

chooseSystem = function(){
  rbVal = as.numeric(tclvalue(rbValue))
  tkdestroy(tt)
  if (rbVal==1)
    tkmessageBox(title = "Sparky Jr. data processing", message="ECD and FID/TCD data.")
  if (rbVal==2)
    tkmessageBox(title = "Sparky III data processing", message="ECD, FID and PDHID data.")
  tclvalue(radiobuttondone) <- 1
  gc <<- rbVal
}# End of function chooseSystem



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function chooseUnits allows user to choose between concentration and
# mass units for GHG amounts.  ppm are preferred for gas samples; ug for
# insitu incubations.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

chooseUnits = function(){
  rbVal = as.numeric(tclvalue(rbValue))
  tkdestroy(tt)
  if (rbVal==1){
    # No longer does anything.  Used to display these cutsie messages...   
    # tkmessageBox(title = "Output units: ppm", message="Good choice, ppm works best for gas samples!")
  }
  if (rbVal==2){
    # tkmessageBox(title = "Output units: ug", message="Good choice, ng GHG is great for insitu incubations!")
  }
  tclvalue(radiobuttondone) <- 1
  units <<- rbVal
}# End of function chooseUnits


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function fitSumTable uses calibration curve to fit area or height data
# according to @var slot from bestCal().  Adds results along with upper
# and lower 95% CI.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


fitSumTable = function(cal.Comp, name_cal, allData){
  print(paste("fitSumTable[", name_cal, "]"))
  nSumCols = ncol(allData)
  insertCol = min(which(names(allData) == paste("Amount_", name_cal, sep = "")), na.rm = TRUE)
  #  comp.name = which(name_cal == fitCompTable[,"cal_name"])
  if (cal.Comp@var == "area"){
    x_data = as.numeric(allData[, paste("Area_", name_cal, sep = "")]) #All Area data
  }else{
    x_data = as.numeric(allData[, paste("Height_", name_cal, sep = "")]) #All Height data
  }
  fit.results = predict(cal.Comp@lm, data.frame(x = x_data), se.fit = TRUE, level = 0.95, interval = "confidence")
  fitTable = cbind(allData[,1:insertCol], fit.results, allData[, (insertCol + 1):nSumCols])
  if (cal.Comp@rsq < 0.95) {
    colnames(fitTable)[insertCol + 1] = paste("bad", colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 1], sep = "_")
    colnames(fitTable)[insertCol + 2] = paste("bad", colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 2], sep = "_")
    colnames(fitTable)[insertCol + 3] = paste("bad", colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 3], sep = "_")
    colnames(fitTable)[insertCol + 4] = paste("bad", colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 4], sep = "_")
    colnames(fitTable)[insertCol + 5] = paste("bad", colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 5], sep = "_")
    colnames(fitTable)[insertCol + 6] = paste("bad", colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 6], sep = "_")
  }else{
    colnames(fitTable)[insertCol + 1] = paste(colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 1], sep = "_")
    colnames(fitTable)[insertCol + 2] = paste(colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 2], sep = "_")
    colnames(fitTable)[insertCol + 3] = paste(colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 3], sep = "_")
    colnames(fitTable)[insertCol + 4] = paste(colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 4], sep = "_")
    colnames(fitTable)[insertCol + 5] = paste(colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 5], sep = "_")
    colnames(fitTable)[insertCol + 6] = paste(colnames(fitTable)[insertCol], colnames(fitTable)[insertCol + 6], sep = "_")
  }
  return(fitTable)
} # End of fitSumTable()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function getCal() generates calibration curve by calling best.lm() and
# interpolates concentrations by calling fitSumTable(). 
#
# Example: getCal("CH4")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# getCal(, SumTable)
# cname = "CO2"
# allData = SumTable
getCal = function(cname, allData){
  cal.lm = new("Regress", lm = lm(1~1), var = "", rsq = 1.0)
  compArea = paste("Area_", cname, sep = "")
  compHeight = paste("Height_", cname, sep = "")
  if ((compArea %in% colnames(allData)) || (compHeight %in% colnames(allData))){ # Check for data columns
    Area_std = as.numeric(allData[which(stds), compArea]) #Stds only
    Height_std = as.numeric(allData[which(stds), compHeight]) #Stds only
    compStd = paste(cname, "_std", sep = "")
    Comp_std = get(compStd)
    Area_std[is.na(Comp_std)] = NA
    Height_std[is.na(Comp_std)] = NA
    Comp_std[is.na(Area_std)] = NA
    Comp_std[is.na(Height_std)] = NA
    stdslist = which(!is.na(Comp_std))
    if (sum(Area_std[!is.na(Comp_std)] > 0, na.rm = TRUE) > 2){ 
      print("Test for minimum of three unique stds in data set")
      gas.conc = Comp_std[stdslist]
      area = Area_std[stdslist]
      height = Height_std[stdslist]
      uni.gas = data.frame(gas.conc, area, height)
      univers = unlist(lapply(c(min.cal:length(gas.conc)), function(x) combn(c(1:length(gas.conc)), x, simplify=FALSE)), recursive=FALSE)
      unitest = unlist(lapply(univers, function(x) length(unique(uni.gas[x, 1])) > min.cal))
      if (sum(unitest=="TRUE") > 0){
        print("use best.lm") #at least one set of three or more unique stds present in data set
        cal.lm = best.lm(compStd, gas.conc, area, height, max.nrmse)
      } else {
        print("use 2-point calibration")
        tkmessageBox(title = paste("2-point calibration for ", cname, sep = ""), message = "Only 2 unique stds present!", icon="warning", type="ok")
        cal.lm = Two_point.lm(compStd, gas.conc, area, height, max.nrmse)
      }
      print(summary(cal.lm@lm))
      print(paste(cname, "calibration residuals:"))
      print(abs(residuals(cal.lm@lm)))
      

    }else print(paste("Insufficient stds to calibrate", cname, "."))
  }else print(paste("No", cname, "data present."))  
  return(cal.lm)
} # End of function getCal()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function getCalTC() generates calibration curve by calling best.lm() and
# interpolates concentrations by calling fitSumTable(). 
# This copy is for TC CO2 data only, based on Area counts.
#
# Example: getCalTC("CH4", SumTable)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 



getCalTC = function(cname, allData){
  cal.lm = new("Regress", lm = lm(1~1), var = "", rsq = 1.0)
  compArea = paste("Area_", cname, sep = "")
  compHeight = paste("Height_", cname, sep = "")
  if (compArea %in% colnames(allData)) { # Check for data columns
    Area_std = as.numeric(allData[which(stds), compArea]) #Stds only
    compStd = paste(cname, "_std", sep = "")
    Comp_std = get(compStd)
    Area_std[is.na(Comp_std)] = NA
    Comp_std[is.na(Area_std)] = NA
    stdslist = which(!is.na(Comp_std))
    if (sum(Area_std[!is.na(Comp_std)] > 0, na.rm = TRUE) > 2){ 
      print("Test for minimum of three unique stds in data set")
      gas.conc = Comp_std[stdslist]
      area = Area_std[stdslist]
      uni.gas = data.frame(gas.conc, area, area)
      univers = unlist(lapply(c(min.cal:length(gas.conc)), function(x) combn(c(1:length(gas.conc)), x, simplify=FALSE)), recursive=FALSE)
      unitest = unlist(lapply(univers, function(x) length(unique(uni.gas[x, 1])) > min.cal))
      if (sum(unitest=="TRUE") > 0){
        print("use best.lm") #at least one set of three or more unique stds present in data set
        cal.lm = best.lm(compStd, gas.conc, area, area, max.nrmse)
      } else {
        print("use 2-point calibration")
        tkmessageBox(title = paste("2-point calibration for ", cname, sep = ""), message = "Only 2 unique stds present!", icon="warning", type="ok")
        cal.lm = Two_point.lm(compStd, gas.conc, area, area, max.nrmse)
      }
      print(summary(cal.lm@lm))
      print(paste(cname, "calibration residuals:"))
      print(abs(residuals(cal.lm@lm)))
      
      
    }else print(paste("Insufficient stds to calibrate", cname, "."))
  }else print(paste("No", cname, "data present."))  
  return(cal.lm)
} # End of function getCalTC()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function lmRegress calculates linear regression model for y versus x,
# plots the graph and returns the lm.y
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

lmRegress = function(y, x, name){
  lm.y = lm(y~x) # Calculate linear regression model
  #corXY = cor.test(x,y) # Pearson r2[4] and p-value[3] for slope
  print(paste("lmRegress: y and x for ", name, ":"))
  print(y)
  print(x)
  print(summary(lm.y))
  x_max = max(x[!is.na(y)], na.rm = TRUE)
  y_max = max(y, na.rm = TRUE)
  return(lm.y)
} # End of function lmRegress()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function plotModel plots linear model with included and excluded data.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


plotModel = function(y, x, model.y, model.name, x.name, y0, x0){
  print(paste("plotModel: [", model.name, "]"))
  x_max = max(x0[!is.na(y0)], na.rm = TRUE)
  y_max = max(y0, na.rm = TRUE)
  corXY = summary(model.y)$r.squared
  windows()
  png(paste(dirRep, paste(model.name, ".png", sep = "" ), sep = ""))
  plot(x0, y0, main = paste("Best fit for", model.name), xlab = x.name, ylab = model.name, xlim = c(0, x_max), ylim = c(0, y_max))
  points(x0, y0, pch = 16, col = "red")
  points(x, y, pch = 16, col = "blue")
  print(y)
  print(x)
  abline(model.y)
  text(0.25*x_max, 0.75*y_max, labels = paste("rsq =", round(as.numeric(corXY), 7)), pos = 4)
  text(0.25*x_max, 0.70*y_max, labels = paste("p-value =",round(as.numeric(cor.test(x, y)[3]), 7)), pos = 4)
  text(0.25*x_max, 0.65*y_max, labels = paste("slope =",round(as.numeric(model.y$coefficients[2]), 7)), pos = 4)
  text(0.25*x_max, 0.60*y_max, labels = paste("inter =",round(as.numeric(model.y$coefficients[1]), 7)), pos = 4)
  dev.off()
  dev.off()
  HTML("<hr>",file=HTMLoutput)
  HTMLInsertGraph(paste(model.name, ".png", sep = ""), file=HTMLoutput, GraphBorder = 3, Align = "center")
  HTML(summary.lm(model.y), file=HTMLoutput)
  HTML("Fitted standard values:", file = HTMLoutput)
  HTML(fitted.values(model.y), file=HTMLoutput)
  HTML("Expected standard values:", file = HTMLoutput)
  HTML((fitted.values(model.y)+residuals(model.y)), file=HTMLoutput)
  if (corXY < 0.95){
    tkmessageBox(title = paste("Low correlation for ", model.name,"!!", sep = ""), message = paste("rsq = ", corXY, sep = ""), icon="warning", type="ok")
  }
} # End of function plotModel()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function setSPKY3Dir() sets up directory assigments for working with
# Sparky 3 data sets.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


setSPKY3Dir = function(){
  dir1 = "F:/Sparky3/SPARKY3SRI/DATA/*.txt"
  dir2 =  "F:/Sparky3/SPARKY3HPCHEM/1/DATA/*.txt"
  testDir1 = "L:/SPARKY3/"
  if (!is.na(file.info(testDir1)$isdir)){
    # Working from T1700
    print(paste("Found", testDir1, ".  Must be working from T1700."))
    dirOut = "L:/SPARKY3/*.csv"
    dirParser = "C:/Martin/data processing"
    dirReport = "L:/SPARKY3/"
    dirWork = "C:/Martin/data processing/R"
  } else {
    testDir1 = "C:/Users/USDA-ARS-Soils526/LabWork/Sparky3/"
    if (!is.na(file.info(testDir1)$isdir)){
      # Working from Lab526?
      print(paste("Found", testDir1, ".  Must be working from Lab526."))
      dir1 = "F:/Rawdata/Sparky3/SPARKY3SRI/DATA/*.txt"
      dir2 =  "F:/Rawdata/Sparky3/SPARKY3HPCHEM/1/DATA/*.txt"
      dirOut = "C:/Users/USDA-ARS-Soils526/LabWork/Sparky3/*.csv"
      dirParser = "C:/Users/USDA-ARS-Soils526/Labwork/Data processing"
      dirReport = "C:/Users/USDA-ARS-Soils526/LabWork/Sparky3/"
      dirWork = "C:/Users/USDA-ARS-Soils526/LabWork/Data processing/R"
    } else {
      # None of the above.  Find your own way.
      print("Could not find any of the expected file paths.  Paths all set to C:/")
      dir1 = "C:/*.txt"
      dir2 = "C:/*.txt"
      dirOut = "C:/*.csv"
      dirWork = "C:/"
      dirParser = "C:/"
    }
  }
  return(new("mydirlist", d1 = dir1, d2 = dir2, dO = dirOut, dP = dirParser, dR = dirReport, dW = dirWork))
}# End of function setSPKY3Dir()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function setSPKYJRDir() sets up directory assigments for working with
# Sparky Jr data sets.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


setSPKYJRDir = function(){
  dir1 = "F:/SparkyJr/SPARKYJRHPCHEM/1/DATA/*.txt"
  dir2 = "F:/SparkyJr/SPARKYJRHPCHEM/2/DATA/*.txt"
  testDir1 = "L:/SPARKYJR/"
  if (!is.na(file.info(testDir1)$isdir)){
    # Working from T1700
    print(paste("Found", testDir1, ".  Must be working from T1700."))
    dirOut = "L:/SPARKYJR/*.csv"
    dirParser = "C:/Martin/data processing"
    dirReport = "L:/SPARKYJR/"
    dirWork = "C:/Martin/data processing/R"
  } else {
    testDir1 = "C:/Users/USDA-ARS-Soils526/LabWork/SparkyJr/"
    if (!is.na(file.info(testDir1)$isdir)){
      # Working from Lab526?
      print(paste("Found", testDir1, ".  Must be working from Lab526."))
      dir1 = "F:/Rawdata/SparkyJr/SPARKYJRHPCHEM/1/DATA/*.txt"
      dir2 = "F:/Rawdata/SparkyJr/SPARKYJRHPCHEM/2/DATA/*.txt"
      dirOut = "C:/Users/USDA-ARS-Soils526/LabWork/SparkyJr/*.csv"
      dirParser = "C:/Users/USDA-ARS-Soils526/Labwork/Data processing"
      dirReport = "C:/Users/USDA-ARS-Soils526/LabWork/SparkyJr/"
      dirWork = "C:/Users/USDA-ARS-Soils526/LabWork/Data processing/R"
    } else {
      # None of the above.  Find your own way.
      print("Could not find any of the expected file paths.  Paths all set to C:/")
      dir1 = "C:/*.txt"
      dir2 = "C:/*.txt"
      dirOut = "C:/*.csv"
      dirWork = "C:/"
      dirParser = "C:/"
    }
  }
  return(new("mydirlist", d1 = dir1, d2 = dir2, dO = dirOut, dP = dirParser, dR = dirReport, dW = dirWork))
}# End of function setSPKYJRDir()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function setSPKYTCDir() sets up directory assigments for working with
# Sparky TC data sets.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


setSPKYTCDir = function(){
  dir1 = ""
  dir2 = "F:/SPARKY/PENData/Data_CO2incubations/*.csv"
  testDir1 = "L:/SPARKY/PENData/Data_CO2incubations/"
  if (!is.na(file.info(testDir1)$isdir)){
    # Working from T1700
    print(paste("Found", testDir1, ".  Must be working from T1700."))
    dirOut = "L:/Sparky/PENData/Data_CO2incubations/*.csv"
    dirParser = "C:/Martin/data processing"
    dirReport = "L:/Sparky/PENData/Data_CO2incubations/"
    dirWork = "C:/Martin/data processing/R"
  } else {
    testDir1 = "C:/Users/USDA-ARS-Soils526/LabWork/SPARKY/PENData/Data_CO2incubations/"
    if (!is.na(file.info(testDir1)$isdir)){
      # Working from Lab526?
      print(paste("Found", testDir1, ".  Must be working from Lab526."))
      dir1 = ""
      dir2 = "F:/Rawdata/SPARKY/PENData/Data_CO2incubations/*.csv"
      dirOut = "C:/Users/USDA-ARS-Soils526/LabWork/Sparky/PENData/Data_CO2incubations/*.csv"
      dirParser = "C:/Users/USDA-ARS-Soils526/Labwork/Data processing"
      dirReport = "C:/Users/USDA-ARS-Soils526/LabWork/Sparky/PENData/Data_CO2incubations/"
      dirWork = "C:/Users/USDA-ARS-Soils526/LabWork/Data processing/R"
    } else {
      # None of the above.  Find your own way.
      print("Could not find any of the expected file paths.  Paths all set to C:/")
      dir1 = ""
      dir2 = "C:/*.csv"
      dirOut = "C:/*.csv"
      dirWork = "C:/"
      dirParser = "C:/"
    }
  }
  return(new("mydirlist", d1 = dir1, d2 = dir2, dO = dirOut, dP = dirParser, dR = dirReport, dW = dirWork))
}# End of function setSPKYTCDir()



# End of function definitions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Select GC system for data to be processed
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

gc = NA
tt <- tktoplevel()
tktitle(tt) <- "Select GC system"
rb1 <- tkradiobutton(tt)
rb2 <- tkradiobutton(tt)
rb3 <- tkradiobutton(tt)
rbValue <- tclVar(1)
tkconfigure(rb1,variable=rbValue,value=1)
tkconfigure(rb2,variable=rbValue,value=2)
tkconfigure(rb3,variable=rbValue,value=3)
tkgrid(tklabel(tt,text="Select GC system:"))
tkgrid(tklabel(tt,text="Sparky Jr."),rb1)
tkgrid(tklabel(tt,text="Sparky III"),rb2)
tkgrid(tklabel(tt,text="Sparky TC"), rb3)
OK.but <- tkbutton(tt, text="OK", command = chooseSystem)
tkgrid(OK.but)
tkfocus(tt)
tkwait.variable(radiobuttondone)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Import Sparky Jr data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

if (gc == 1){
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Sparky Jr. specific constants
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  systemName = "SPKYJR"
  TCDList = "CO2N2"
  
  mydir = setSPKYJRDir()
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Select ECD and FID data files to be parsed.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  #ECD data file contains sample list and must be processed first.
  cat("Imported ECD data from:\n\n")
  
  infile2 = choose.files(mydir@d2, filters = myFilters[c("txt","All"),], caption = "Choose SPARKY JR. ECD datafile")
  
  ECDData = read.csv(infile2, sep = " ", header = FALSE, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  cat(infile2, "\n\n")
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Setup output directory based on sequence name.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  sequenceName = strsplit(infile2, "\\\\")[[1]][length(strsplit(infile2, "\\\\")[[1]]) - 1]
  dirRep = paste(mydir@dR, sequenceName, "/", sep="")
  if(is.na(file.info(dirRep)$isdir) == FALSE){
    #Output folder exists.  Do nothing
  }else{
    #Output folder does not exist.  Create new folder.
    dir.create(dirRep)
  }
  
  HTMLoutput=file.path(dirRep,"pReport.html")
  HTML.title(paste("Regression summary data for ", sequenceName, " (", systemName, ")", sep=""), Align = "center", HR=3, file=HTMLoutput)
  HTML.title(Sys.time(), Align = "center", HR=4, file=HTMLoutput)
  
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Locate first Compound Summary Table in data file and use the column
  # separator character "|" to identify column widths.  Column widths
  # are stored in Seqfwf and used by read.fwf()to extract sequence
  # information: Run, Vial, Inj, Date and Time, Filename.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  runECD = which(ECDData[, 1] == "Run")
  fw = runECD + 2
  l = ECDData[fw[1],1] 
  m=strsplit(l,"|")
  Seqfw = which(m[[1]][] == "|")
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Manual fix for start of File Name column
  
  Seqfw[4] = Seqfw[4] - 1
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Seqfwf = Seqfw[1]
  for (i in 1:length(Seqfw)-1){
    newSeqfw = Seqfw[i+1] - Seqfw[i]
    Seqfwf = c(Seqfwf, newSeqfw)
  }
  lastSeqfw = length(m[[1]]) - Seqfw[i+1]
  Seqfwf = c(Seqfwf, lastSeqfw)
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Now re-read ECD data with correct fwf for sequence information
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  Sequence = read.fwf(infile2, widths = Seqfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  runSeq = which(Sequence[, 1] == "Run")
  nSamples = runSeq[2] - runSeq[1] - 6
  nSeqCols = ncol(Sequence)
  beginSeq = runSeq[1] + 2
  endSeq = runSeq[2] - 5
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Determine fwf info for compound information in ECD datafile:
  # Run, Type, RetTime, Amount, Area, Height, Width, Symm.
  # Store in ECDfwf.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  l = ECDData[fw[2],1] 
  m=strsplit(l,"|")
  ECDfw = as.numeric(which(m[[1]][] == "|"))
  ECDfwf = ECDfw[1]
  for (i in 1:length(ECDfw)-1){
    newECDfw = ECDfw[i+1]-ECDfw[i]
    ECDfwf = c(ECDfwf, newECDfw)
  }
  lastECDfw = length(m[[1]]) - ECDfw[i+1]
  ECDfwf = c(ECDfwf, lastECDfw)
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Re-read ECD datafile using ECDfwf to extract compound information.
  # Overwrite ECDData array with fwf compound data.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ECDData = read.fwf(infile2, widths = ECDfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  runECD = which(ECDData[, 1] == "Run")
  endECD = which(ECDData[, 1] == "Mean")
  nECDCols = ncol(ECDData)
  nCompECD = length(endECD)
  beginECD = runECD + 2
  endECD = endECD - 2
  compECD = runECD - 2
  
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Now get fwf for FID data
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  # Try to locate FIDData automatically
  print("Try to find FIDData file automatically.")
  infile1 = infile2
  if ("\\\\2" %in% infile1) { # ECDData file located in HPCHEM folder structure
  infile1 = sub("\\\\2", "\\\\1", infile1)
  } else { # ECDData file not located in HPCHEM folder structure
  infile1 = sub("\\ECDFID.txt", "\\FIDTCD.txt", infile1)
  }
  if (file.exists(infile1)){
    cat("Importing FIDTCD data.")
  } else {
    infile1 = choose.files(mydir@d1, filters = Filters[c("txt","All"),], caption = "Choose SPARKY JR. FID datafile")
  }
  FIDData = read.csv(infile1, sep = " ", header = FALSE, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  cat("Imported FID data from:\n\n")
  cat(infile1, "\n\n")
  runFID = which(FIDData[, 1] == "Run")
  fw = runFID + 2
  l = FIDData[fw[2],1] 
  m=strsplit(l,"|")
  FIDfw = as.numeric(which(m[[1]][] == "|")
                     )
  FIDfwf = FIDfw[1]
  for (i in 1:length(FIDfw)-1){
    newFIDfw = FIDfw[i+1]-FIDfw[i]
    FIDfwf = c(FIDfwf, newFIDfw)
  }
  lastFIDfw = length(m[[1]]) - FIDfw[i+1]
  FIDfwf = c(FIDfwf, lastFIDfw)
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Now re-read FID data with correct fwf for compound information
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  FIDData = read.fwf(infile1, widths = FIDfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  runFID = which(FIDData[, 1] == "Run")
  endFID = which(FIDData[, 1] == "Mean")
  nFIDCols = ncol(FIDData)
  nCompFID = length(endFID)
  beginFID = runFID + 2
  endFID = endFID - 2
  compFID = runFID - 2
  nSumCols = nSeqCols + nECDCols * nCompECD + nFIDCols * nCompFID + 1
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  dataHeader = matrix(nrow = 1, ncol = nSumCols)
  SumTable = matrix(nrow = nSamples, ncol = nSumCols)
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  headerRow = min(runSeq)
  dataHeader = Sequence[headerRow, ]
  SumTable = Sequence[beginSeq:endSeq, ]
  
  for (k in 1:nCompECD) {
    startCol = ncol(SumTable)
    for (x in beginECD[k+1]:endECD[k]) {
      i = as.integer(ECDData[x, 1])
      for (j in 1:nECDCols) {
        SumTable[i, (j+startCol)] = ECDData[x, j]
      }
    }
    for (j in 1:nECDCols)
      dataHeader[(j+startCol)] =  paste(ECDData[runECD[k+1],j], strsplit(ECDData[compECD[k+1], 3], " ")[[1]][1], "ECD", sep = "_")
  }
  for (k in 1:nCompFID) {
    startCol = ncol(SumTable)
    for (x in beginFID[k+1]:endFID[k]) {
      i = as.integer(FIDData[x, 1])
      for (j in 1:nFIDCols) {
        SumTable[i, (j+startCol)] = FIDData[x, j]
      }
    }
    for (j in 1:nFIDCols) {
      compName = strsplit(FIDData[compFID[k+1], 3], " ")[[1]][1]
      if (length(grep(compName, TCDList, ignore.case = TRUE)) != 0)
        dataHeader[(j+startCol)] =  paste(FIDData[runFID[k+1],j], compName, sep = "_")
      else
        dataHeader[(j+startCol)] =  paste(FIDData[runFID[k+1],j], compName, sep = "_")
    }
  }
  
  colnames(SumTable) = dataHeader
  in1 = paste("...", strsplit(toupper(infile1), "DATA")[[1]][2], sep = "")
  in2 = paste("...", strsplit(toupper(infile2), "DATA")[[1]][2], sep = "")
  tkmessageBox(title="SparkyJr: Parsing successful!", message=paste("Parsed files:", in2,"and", in1), icon="info", type="ok")
} # End of SparkyJr data import

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Import Sparky III data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

if (gc == 2){
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Sparky III specific constants
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  npd = 3  # Number of generic columns in PD HID data file
  cName = 4 # Column number of first compound name
  pdCols = 7	# Number of columns per compound in PD HID data file
  systemName = "SPKY3"
  ECDList = "CO2N2O"
  
  mydir = setSPKY3Dir()
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Select ECD/FID data file to be parsed.
  # ECD/FID data file contains sample list and must be processed first.
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  cat("Imported ECD/FID data from:\n\n")
  
  infile2 = choose.files(mydir@d2, filters = myFilters[c("txt","All"),], caption = "Choose SPARKY 3 ECD/FID datafile")
  ECDData = read.csv(infile2, sep = " ", header = FALSE, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  #edit(ECDData)
  cat(infile2, "\n\n")
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Setup output directory based on sequence name.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  sequenceName = strsplit(infile2, "\\\\")[[1]][length(strsplit(infile2, "\\\\")[[1]]) - 1]
  dirRep = paste(mydir@dR, sequenceName, "/", sep="")
  if(is.na(file.info(dirRep)$isdir) == FALSE){
    #Output folder exists.  Do nothing
  }else{
    #Output folder does not exist.  Create new folder.
    dir.create(dirRep)
  }
  
  HTMLoutput=file.path(dirRep,"pReport.html")
  HTML.title(paste("Regression summary data for ", sequenceName, " (", systemName, ")", sep=""), Align = "center", HR=3, file=HTMLoutput)
  HTML.title(Sys.time(), Align = "center", HR=4, file=HTMLoutput)
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Locate first Compound Summary Table in data file and use the column
  # separator character "|" to identify column widths.  Column widths
  # are stored in Seqfwf and used by read.fwf()to extract sequence
  # information: Run, Vial, Inj, Date and Time, Filename.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  runECD = which(ECDData[, 1] == "Run")
  fw = runECD + 2
  l = ECDData[fw[1],1] 
  m=strsplit(l,"|")
  Seqfw = which(m[[1]][] == "|")
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Manual fix for start of File Name column
  
  Seqfw[4] = Seqfw[4] - 1
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Seqfwf = Seqfw[1]
  for (i in 1:length(Seqfw)-1){
    newSeqfw = Seqfw[i+1]-Seqfw[i]
    Seqfwf = c(Seqfwf, newSeqfw)
  }
  lastSeqfw = length(m[[1]]) - Seqfw[i+1]
  Seqfwf = c(Seqfwf, lastSeqfw)  #Sequence fixed width information
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Re-read ECD datafile using Seqfwf to extract sequence information.
  # Store sequence information in Sequence array.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Sequence = read.fwf(infile2, widths = Seqfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  runSeq = which(Sequence[, 1] == "Run")
  nSamples = runSeq[2] - runSeq[1] - 6
  nSeqCols = ncol(Sequence)
  beginSeq = runSeq[1] + 2
  endSeq = runSeq[2] - 5
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Determine fwf info for compound information in ECD datafile:
  # Run, Type, RetTime, Amount, Area, Height, Width, Symm.
  # Store in ECDfwf.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  l = ECDData[fw[2],1] 
  m=strsplit(l,"|")
  ECDfw = as.numeric(which(m[[1]][] == "|"))
  ECDfwf = ECDfw[1]
  for (i in 1:length(ECDfw)-1){
    newECDfw = ECDfw[i+1]-ECDfw[i]
    ECDfwf = c(ECDfwf, newECDfw)
  }
  lastECDfw = length(m[[1]]) - ECDfw[i+1]
  ECDfwf = c(ECDfwf, lastECDfw)
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Re-read ECD datafile using ECDfwf to extract compound information.
  # Overwrite ECDData array with fwf compound data.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ECDData = read.fwf(infile2, widths = ECDfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Locate beginning and end of each compound summary table in ECDData array.
  # Count the number of different compounds identified in the datafile.
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  runECD = which(ECDData[, 1] == "Run")
  endECD = which(ECDData[, 1] == "Mean")
  runVal = which(is.na(ECDData[, 1]) == "TRUE")
  maxRunPos = runVal[min(which(runVal > runECD[1]))] - 1
  nSamples = as.numeric(ECDData[maxRunPos, 1])
  nECDCols = ncol(ECDData)  #number of columns in ECDData array
  nCompECD = length(endECD) #number of compounds summarized in ECDData array; same as number of compuounds in data file
  beginECD = runECD + 2     #array of first sample rows for each compound summary table
  endECD = endECD - 2       #array of last sample rows for each compound summary table
  compECD = runECD - 2      #array of rows containing compound name for each compound summary table
  nSumCols = nSeqCols + nECDCols * nCompECD + 1  #total number of columns in Summary Table: SumTable
  
  dataHeader = matrix(nrow = 1, ncol = nSumCols)
  SumTable = matrix(nrow = nSamples, ncol = nSumCols)
  headerRow = min(runSeq)
  dataHeader = Sequence[headerRow, ]
  SumTable = Sequence[beginSeq:maxRunPos, ]
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Fill Summary Table array using Run number as row index i.
  # Generate column names for SumTable.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  for (k in 1:nCompECD) {
    startCol = ncol(SumTable)
    for (z in beginECD[k+1]:endECD[k]) {
      i = as.integer(ECDData[z, 1])
      for (j in 1:nECDCols) {
        SumTable[i, (j+startCol)] = ECDData[z, j]
      }
    }
    for (j in 1:nECDCols) {
      dataHeader[(j+startCol)] = ECDData[runECD[k+1], j]
      if (ECDData[runECD[k+1], j] == "Amount") {
        compInfo = ""
        for (l in 1:3){  #Only paste the first 3 columns
          compInfo = paste(compInfo, ECDData[compECD[k+1], l], sep = "")
        }
      }
    }
    for (j in 1:nECDCols) {
      compName = strsplit(ECDData[compECD[k+1], 3], " ")[[1]][1]
      if (length(grep(compName, ECDList, ignore.case = TRUE)) != 0)
        dataHeader[(j+startCol)] =  paste(ECDData[runECD[k+1],j], compName, "ECD", sep = "_")
      else
        # CH4 is the other compound in this data file.  Do not add detector identifier
        dataHeader[(j+startCol)] =  paste(ECDData[runECD[k+1],j], compName, sep = "_")
    }
  }
  
  colnames(SumTable)= dataHeader
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Now get fwf for PDHID data
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  cat("Imported SRI data from:\n\n")
  infile1 = mydir@d1
  infile1 = sub("\\*", sequenceName, infile1)
  if (file.exists(infile1)){
    cat("Importing SRI data.")
  } else {
    infile1 = choose.files(mydir@d1, filters = myFilters[c("txt","All"),], caption = "Choose SPARKY3 PDHID/SRI datafile")
  }
  nf = max(count.fields(infile1, sep = "\t")) 
  pdComp = (nf - npd)/pdCols
  PDHIDData = read.csv(infile1, sep = "\t", header = FALSE, fill = TRUE)
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  # Manage PDHIDData headers
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  pdHeader = matrix(nrow = 1, ncol = nf)
  pdHeader[1,1:npd] = c("Filename_PID", "Date_PID", "Time_PID")
  npd = npd + 1
  for (p in 1:pdComp - 1) {
    mpd = npd + pdCols - 1
    pdHeader[1,npd:mpd] = paste(c("Compound", "RT", "Amount", "Area", "Height", "Width", "Symmetry"), PDHIDData[1, cName], sep  = "_")
    cName = cName + pdCols
    npd = mpd + 1
  }
  
  # Append PID identifier to CH4 and N2O data columns to distinguish from
  # FID and ECD data columns
  
  CH4_PDHID = grep("CH4", pdHeader, ignore.case = TRUE)
  pdHeader[CH4_PDHID] = paste(pdHeader[CH4_PDHID], "PID", sep = "_") 
  N2O_PDHID = grep("N2O", pdHeader, ignore.case = TRUE)
  pdHeader[N2O_PDHID] = paste(pdHeader[N2O_PDHID], "PID", sep = "_") 
  
  if(ncol(pdHeader) == ncol(PDHIDData)){
    colnames(PDHIDData) = pdHeader
  }else {
    tkmessageBox(title="Sparky3: Parser Error 103", message=paste("Error 103: PDHID/SRI Column number mismatch in", basename(infile1),". Look for double or missing peak in data file."), icon="warning", type="ok")
    parserError = "Error 103"
  }
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Column bind SumTable and PDHIDData into final summary array: SumTable
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  in1 = paste("...", strsplit(toupper(infile1), "\\DATA")[[1]][2], sep = "")
  in2 = paste("...", strsplit(toupper(infile2), "\\DATA")[[1]][2], sep = "")
  
  if (nrow(SumTable)  == nrow(PDHIDData)){
    SumTable = cbind(SumTable, PDHIDData)
    tkmessageBox(title="Sparky3: Parsing successful!", message=paste("Parsed files:", in2,"and", basename(infile1)), icon="info", type="ok")
  }else{
    tkmessageBox(title="Sparky3: Parser Error 123", message=paste("Error 123: Row number mismatch between", in2," [",nrow(SumTable),"] and", basename(infile1)," [",nrow(PDHIDData),"]"), icon="warning", type="ok")
    parserError = "Error 123"
  }
  
  cat("parser Error:",parserError, "\n")
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  # End of Sparky III data import
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Import Sparky TotalChrom data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

if (gc == 3){
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Sparky TotalChrom specific constants
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  systemName = "SPKYTC"
  compName = ""
  TCheader = c(NA)
  DateTime = c(NA)
  
  mydir = setSPKYTCDir()
 
  infile2 = choose.files(mydir@d2, filters = myFilters, caption = "Import SPARKY FID/TCD data")
  TCData = read.csv(infile2, sep = ",", header = FALSE, as.is = TRUE, strip.white = TRUE)
  cat(infile2, "\n\n")  
  
  sequenceName = strsplit(infile2, "\\\\")[[1]][length(strsplit(infile2, "\\\\")[[1]]) - 1]
  dirRep = paste(mydir@dR, sequenceName, "/", sep="")
  if(is.na(file.info(dirRep)$isdir) == FALSE){
    #Output folder exists.  Do nothing
  }else{
    #Output folder does not exist.  Create new folder.
    dir.create(dirRep)
  }
  
  HTMLoutput=file.path(dirRep,"pReport.html")
  HTML.title(paste("Regression summary data for ", sequenceName, " (", systemName, ")", sep=""), Align = "center", HR=3, file=HTMLoutput)
  HTML.title(Sys.time(), Align = "center", HR=4, file=HTMLoutput)
  

  for (i in 1:(ncol(TCData)-1)){
    str1 = min(which(TCData[,i]!=""), na.rm = TRUE)
    if (str1 == 2) {
      compName = ""
      TCheader[i] = paste(compName, TCData[str1,i],TCData[str1+1,i])
      compName = substr(TCData[str1,i], 1, 3)
    }else if ("Area" %in% TCData[str1,i]){
      TCheader[i] = paste(TCData[str1,i], compName, sep = "_")
    }else if("Amount" %in% TCData[str1+1,i]){
      TCheader[i] = paste(TCData[str1+1,i], compName, sep = "_")
    }else{
      TCheader[i] = paste(compName, TCData[str1,i],TCData[str1+1,i])
    }
  }
  
  
  
  TCheader = sub("^\\s+", "", TCheader)
  colnames(TCData) = TCheader
  
  # colnames(TCData) = sub("^\\s+", "", TCheader) # Remove leading blanks from column names
  
  
  # sub("^\\s+", "", " V 3 d")
  # "File Name"                            
  # [3] "Date of Injection"                     "Time of Injection"                    
  # [5] "Sample Name"                           "Sample Number"                        
  # [7] "AIR TCD Time"                          "Area [µV·s]"                          
  # [9] "Adjusted Amount"                       "CH4 TCD Time"                         
  # [11] "Area [µV·s]"                           "Adjusted Amount"                      
  # [13] "CO2 TCD ug-C Time"                     "Area [µV·s]"                          
  # [15] "Adjusted Amount"                       "NA NA"     
  
  
  # Delete all rows prior to sample data.
  
  # Maybe use this to only export rows with data to summary file
  
  startRow = min(grep(".rst", TCData[, "File Name"], ignore.case = TRUE), na.rm = TRUE)
  endRow = max(grep(".rst", TCData[, "File Name"], ignore.case = TRUE), na.rm = TRUE)
  SumTable = TCData[startRow:endRow,]
  nSamples = length(grep(".rst", SumTable[, "File Name"], ignore.case = TRUE))
  colnames(SumTable)
  DoI = grep("Date of Injection",colnames(SumTable), ignore.case = FALSE)
  ToI = grep("Time of Injection",colnames(SumTable), ignore.case = FALSE)
  listCols = which(1:ncol(SumTable) != c(DoI, ToI))
  DateTime = paste(SumTable[,"Date of Injection"], SumTable[,"Time of Injection"])
  SumTable = cbind(SumTable[,1:2],DateTime, SumTable[,3:ncol(SumTable)])
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  # Now write SumTable to .csv file
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
#   
#   outfile = choose.files(mydir@dO, filters = myFilters, caption = "Choose SPKYNetFit output file")
#   write.csv(SumTable, file = outfile, quote = FALSE, row.names = FALSE)
#   
  
}



  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  # Choose whether to perform automated or manual calibration
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  calmode = NA
  tt <- tktoplevel()
  tktitle(tt) <- "Select Auto or Manual calibration"
  rb1 <- tkradiobutton(tt)
  rb2 <- tkradiobutton(tt)
  rbValue <- tclVar(1)
  tkconfigure(rb1,variable=rbValue,value=1)
  tkconfigure(rb2,variable=rbValue,value=2)
  tkgrid(tklabel(tt,text="Select calibration mode:"))
  tkgrid(tklabel(tt,text="Auto"),rb1)
  tkgrid(tklabel(tt,text="Manual"),rb2)
  OK.but <- tkbutton(tt, text="OK", command = chooseCalMode)
  tkgrid(OK.but)
  tkfocus(tt)
  tkwait.variable(radiobuttondone)
  
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #  Perform auto-calibration of gc data.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  if (calmode == 1) {
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # If parsing was successful (parserError = "None") then proceed to calculate
    # least squares linear fit of std area and height, select best fit and generate
    # sample values based best fit.
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if (parserError == "None"){
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Rename methane, acetylene, ethylene, ethane column headers to contain
      # chemical formula (CH4) instead of the name.
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      MethCol = grep("Methane",colnames(SumTable), ignore.case = TRUE) # All FID data prior to 091813 had 2 sets of methane data
      if (length(MethCol) > 0){
        for (i in 1:length(MethCol)){
          namesplit = strsplit(colnames(SumTable)[MethCol[i]], "_")
          colnames(SumTable)[MethCol[i]] = paste(namesplit[[1]][1], "CH4", sep = "_")
        }
      }
      
      MethCol = grep("Acety",colnames(SumTable), ignore.case = TRUE)
      if (length(MethCol) > 0){
        for (i in 1:length(MethCol)){
          namesplit = strsplit(colnames(SumTable)[MethCol[i]], "_")
          colnames(SumTable)[MethCol[i]] = paste(namesplit[[1]][1], "C2H2", sep = "_")
        }
      }
      
      MethCol = grep("Ethyl",colnames(SumTable), ignore.case = TRUE)
      if (length(MethCol) > 0){
        for (i in 1:length(MethCol)){
          namesplit = strsplit(colnames(SumTable)[MethCol[i]], "_")
          colnames(SumTable)[MethCol[i]] = paste(namesplit[[1]][1], "C2H4", sep = "_")
        }
      }
      
      MethCol = grep("_Etha",colnames(SumTable), ignore.case = TRUE)
      if (length(MethCol) > 0){
        for (i in 1:length(MethCol)){
          namesplit = strsplit(colnames(SumTable)[MethCol[i]], "_")
          colnames(SumTable)[MethCol[i]] = paste(namesplit[[1]][1], "C2H6", sep = "_")
        }
      }
      
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Compound calibration and regression: Use all available stds
      # to generate the best fit based on 3 (min.cal) or more values.
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Extract stds based on name fragment and populate std data arrays for each compound
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      samplenames = SumTable[,"Sample Name"]
      stdnames = c("HE","5K", "3K", "AMB", "1%", "10%", "AC", "AIR")
      qaqcnames = c("T", "TR", "TS", "Q", "QA", "QC", "-AMB")
      
      # stds contains list of stds in datafile
      
      stds = matrix(nrow = nSamples, ncol = 1)
      for (i in 1:nSamples){
        for (j in 1:length(stdnames)){
          if (length(grep(stdnames[j], SumTable[,"Sample Name"][i], ignore.case = TRUE)) > 0)
            stds[i] = TRUE
        }
      }
      
      # remove possible trip standards or qa/qc standards from list
      
      for (i in 1:nSamples){
        for (j in 1:length(qaqcnames)){
          if (length(grep(qaqcnames[j], SumTable[,"Sample Name"][i], ignore.case = TRUE)) > 0)
            stds[i] = FALSE
        }
      }
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Select concentration or mass units for calibration and calculations using radio button
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      radiobuttondone <- tclVar(0)  
      units = NA
      tt <- tktoplevel()
      tktitle(tt) <- "Select GHG units"
      rb1 <- tkradiobutton(tt)
      rb2 <- tkradiobutton(tt)
      rbValue <- tclVar(1)
      tkconfigure(rb1,variable=rbValue,value=1)
      tkconfigure(rb2,variable=rbValue,value=2)
      tkgrid(tklabel(tt,text="Which units do you prefer?"))
      tkgrid(tklabel(tt,text="ppm & ppb"),rb1)
      tkgrid(tklabel(tt,text=" ug & ng "),rb2)
      OK.but <- tkbutton(tt,text="OK", command = chooseUnits)
      tkgrid(OK.but)
      tkfocus(tt)
      tkwait.variable(radiobuttondone)  
      
      
      if (units == 1) {
        HE = grep(stdnames[1],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[HE] = NA
        N2O_PID_std[HE] = NA
        N2O_ECD_std[HE] = NA
        CH4_std[HE] = NA
        CH4_PID_std[HE] = NA
        C2H2n_std[HE] = NA
        CO2_ECD_std[HE] = NA
        O2_std[HE] = NA
        N2_std[HE] = NA
        
        K5 = grep(stdnames[2],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[K5] = 5000.0
        N2O_PID_std[K5] = 5000.0
        N2O_ECD_std[K5] = 5000.0    
        CH4_std[K5] = 5.0
        CH4_PID_std[K5] = 5.0
        C2H2n_std[K5] = 2.5
        CO2_ECD_std[K5] = 0.5
        O2_std[K5] = 20.0
        N2_std[K5] = 80.0
        
        K3 = grep(stdnames[3],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[K3] = 3000.0
        N2O_PID_std[K3] = 800.0
        N2O_ECD_std[K3] = 800.0
        CH4_std[K3] = 2.0
        CH4_PID_std[K3] = 2.0
        C2H2n_std[K3] = 1.0
        CO2_ECD_std[K3] = 0.3
        O2_std[K3] = 10.0
        N2_std[K3] = 90.0
        
        AMB = grep(stdnames[4],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[AMB] = 400.0
        N2O_PID_std[AMB] = 200.0
        N2O_ECD_std[AMB] = 200.0
        CH4_std[AMB] = 1.0
        CH4_PID_std[AMB] = 1.0
        C2H2n_std[AMB] = 0.5
        CO2_ECD_std[AMB] = 0.04
        O2_std[AMB] = 5.0
        N2_std[AMB] = NA
        
        PCNT1 = grep(stdnames[5],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[PCNT1] = 10000.0
        N2O_PID_std[PCNT1] = NA
        N2O_ECD_std[PCNT1] = NA
        CH4_std[PCNT1] = NA
        CH4_PID_std[PCNT1] = NA
        C2H2n_std[PCNT1] = NA
        CO2_ECD_std[PCNT1] = 1.0
        O2_std[PCNT1] = 2.0
        N2_std[PCNT1] = 1.0
        
        PCNT10 = grep(stdnames[6],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[PCNT10] = NA
        N2O_PID_std[PCNT10] = NA
        N2O_ECD_std[PCNT10] = NA
        CH4_std[PCNT10] = NA
        CH4_PID_std[PCNT10] = NA
        C2H2n_std[PCNT10] = NA
        CO2_ECD_std[PCNT10] = 10.0
        O2_std[PCNT10] = NA
        N2_std[PCNT10] = NA
        
        AC = grep(stdnames[7],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[AC] = NA
        N2O_PID_std[AC] = NA
        N2O_ECD_std[AC] = NA
        CH4_std[AC] = NA
        CH4_PID_std[AC] = NA
        C2H2n_std[AC] = NA
        CO2_ECD_std[AC] = NA
        O2_std[AC] = NA
        N2_std[AC] = NA
        
        AIR = grep(stdnames[8],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[AIR] = NA
        N2O_PID_std[AIR] = NA
        N2O_ECD_std[AIR] = NA
        CH4_std[AIR] = NA
        CH4_PID_std[AIR] = NA
        C2H2n_std[AIR] = NA
        CO2_ECD_std[AIR] = NA
        O2_std[AIR] = NA
        N2_std[AIR] = NA
      }else if (units == 2){
        HE = grep(stdnames[1],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[HE] = NA
        N2O_PID_std[HE] = NA
        N2O_ECD_std[HE] = NA
        CH4_std[HE] = NA
        CH4_PID_std[HE] = NA
        C2H2n_std[HE] = NA
        CO2_ECD_std[HE] = NA
        O2_std[HE] = NA
        N2_std[HE] = NA
        
        K5 = grep(stdnames[2],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[K5] = 12.27
        N2O_PID_std[K5] = 28.63
        N2O_ECD_std[K5] = 28.63
        CH4_std[K5] = 12.27
        CH4_PID_std[K5] = 12.27
        C2H2n_std[K5] = 12.27
        CO2_ECD_std[K5] = 12.27
        O2_std[K5] = 20.0
        N2_std[K5] = 80.0
        
        K3 = grep(stdnames[3],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[K3] = 7.36
        N2O_PID_std[K3] = 4.58
        N2O_ECD_std[K3] = 4.58
        CH4_std[K3] = 4.91
        CH4_PID_std[K3] = 4.91
        C2H2n_std[K3] = 4.91
        CO2_ECD_std[K3] = 7.36
        O2_std[K3] = 10.0
        N2_std[K3] = 90.0
        
        AMB = grep(stdnames[4],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[AMB] = 0.98
        N2O_PID_std[AMB] = 1.15
        N2O_ECD_std[AMB] = 1.15
        CH4_std[AMB] = 2.45
        CH4_PID_std[AMB] = 2.45
        C2H2n_std[AMB] = 2.45
        CO2_ECD_std[AMB] = 0.98
        O2_std[AMB] = 5.0
        N2_std[AMB] = NA
        
        PCNT1 = grep(stdnames[5],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[PCNT1] = 24.54
        N2O_PID_std[PCNT1] = NA
        N2O_ECD_std[PCNT1] = NA
        CH4_std[PCNT1] = NA
        CH4_PID_std[PCNT1] = NA
        C2H2n_std[PCNT1] = NA
        CO2_ECD_std[PCNT1] = 24.54
        O2_std[PCNT1] = 2.0
        N2_std[PCNT1] = 1.0	
        
        PCNT10 = grep(stdnames[6],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[PCNT10] = NA
        # Normally the 10% std is excluded from the calibration because it skews the fit at lower
        # concentrations.  It is included in TC datasets because those usually have very high CO2 concentrations
        if (gc == 3) CO2_std[PCNT10] = 245.4 
        N2O_PID_std[PCNT10] = NA
        N2O_ECD_std[PCNT10] = NA
        CH4_std[PCNT10] = NA
        CH4_PID_std[PCNT10] = NA
        C2H2n_std[PCNT10] = NA
        CO2_ECD_std[PCNT10] = 245.4
        O2_std[PCNT10] = NA
        N2_std[PCNT10] = NA
        
        AC = grep(stdnames[7],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[AC] = NA
        N2O_PID_std[AC] = NA
        N2O_ECD_std[AC] = NA
        CH4_std[AC] = NA
        CH4_PID_std[AC] = NA
        C2H2n_std[AC] = NA
        CO2_ECD_std[AC] = NA
        O2_std[AC] = NA
        N2_std[AC] = NA
        
        AIR = grep(stdnames[8],SumTable[which(stds), "Sample Name"], ignore.case = TRUE)
        CO2_std[AIR] = NA
        N2O_PID_std[AIR] = NA
        N2O_ECD_std[AIR] = NA
        CH4_std[AIR] = NA
        CH4_PID_std[AIR] = NA
        C2H2n_std[AIR] = NA
        CO2_ECD_std[AIR] = NA
        O2_std[AIR] = NA
        N2_std[AIR] = NA
      }
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # 1. Check for data
      # 2. Calculate calibrations for compounds
      # 3. Generate predicted values with confidence interval and insert into SumTable
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # CH4
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      ch4cal = getCal("CH4", SumTable)
      if(ch4cal@var != "") SumTable = fitSumTable(ch4cal, "CH4", SumTable)
      
      ch4pidcal = getCal("CH4_PID", SumTable)
      if(ch4pidcal@var != "") SumTable = fitSumTable(ch4pidcal, "CH4_PID", SumTable)
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Acetylene, ethane, ethylene
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      # Use the CH4 calibration as the basis of C2H2n calibrations
      
      if (("Area_CH4" %in% colnames(SumTable)) || ("Height_CH4" %in% colnames(SumTable))){ # Check for data columns
        Area_C2H2n_std = as.numeric(SumTable[which(stds), "Area_CH4"])*0.5 #Stds only
        Height_C2H2n_std = as.numeric(SumTable[which(stds), "Height_CH4"])*0.5 #Stds only
        C2H2n_std = CH4_std
        Area_C2H2n_std[is.na(C2H2n_std)] = NA
        Height_C2H2n_std[is.na(C2H2n_std)] = NA
        C2H2n_std[is.na(Area_C2H2n_std)] = NA
        C2H2n_std[is.na(Height_C2H2n_std)] = NA
        stdslist = which(!is.na(C2H2n_std))
        if (sum(Area_C2H2n_std[!is.na(C2H2n_std)] > 0, na.rm = TRUE) > 2){ 
          gas.name = "C2H2n_std"
          gas.conc = C2H2n_std[stdslist]
          area = Area_C2H2n_std[stdslist]
          height = Height_C2H2n_std[stdslist]
          cal.C2H2n = best.lm(gas.name, gas.conc, area, height, max.nrmse)
          print(summary(cal.C2H2n@lm))
          print("C2H2n calibration residuals:")
          print(abs(residuals(cal.C2H2n@lm)))
        }else{print("Insufficient stds to calibrate C2H2n!")}
      }else{print("No C2H2n data present.")}
      
      # Use the C2H2n calibration for all the C2-compounds
      
      if (exists("cal.C2H2n")){
        cal.C2H2 = cal.C2H2n 
        cal.C2H4 = cal.C2H2n 
        cal.C2H6 = cal.C2H2n 
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        # C2H2
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        
        if (("Area_C2H2" %in% colnames(SumTable)) || ("Height_C2H2" %in% colnames(SumTable))){
          nSumCols = ncol(SumTable)
          insertCol = which(names(SumTable) == "Amount_C2H2") # Locate column by column name
          SumTable = fitSumTable(cal.C2H2, "C2H2", SumTable)
        }else print("No C2H2 data present.")
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        # C2H4
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        if (("Area_C2H4" %in% colnames(SumTable)) || ("Height_C2H4" %in% colnames(SumTable))){    
          nSumCols = ncol(SumTable)
          insertCol = which(names(SumTable) == "Amount_C2H4") # Locate column by column name
          SumTable = fitSumTable(cal.C2H4, "C2H4", SumTable)
        }else print("No C2H4 data present.")
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        # C2H6
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        
        if (("Area_C2H6" %in% colnames(SumTable)) || ("Height_C2H6" %in% colnames(SumTable))){    
          nSumCols = ncol(SumTable)
          insertCol = which(names(SumTable) == "Amount_C2H6") # Locate column by column name
          SumTable = fitSumTable(cal.C2H6, "C2H6", SumTable)
        }else print("No C2H6 data present.")
      }else print("Insufficient stds to calibrate C2H2n compounds!")
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Generate calibration curves and interpolate data for N2O, N2O_ECD, CO2,
      # CO2_ECD, O2, and N2.
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      n2oecdcal = getCal("N2O_ECD", SumTable)
      if(n2oecdcal@var != "") SumTable = fitSumTable(n2oecdcal, "N2O_ECD", SumTable)
      n2opidcal = getCal("N2O_PID", SumTable)
      if(n2opidcal@var != "") SumTable = fitSumTable(n2opidcal, "N2O_PID", SumTable)
      co2ecdcal = getCal("CO2_ECD", SumTable)
      if(co2ecdcal@var != "") SumTable = fitSumTable(co2ecdcal, "CO2_ECD", SumTable)
      if (gc == 3){        
        co2cal = getCalTC("CO2", SumTable) #Area-only version of getCal
      } else {
        co2cal = getCal("CO2", SumTable)
      }
      if(co2cal@var != "") SumTable = fitSumTable(co2cal, "CO2", SumTable)
      o2cal = getCal("O2", SumTable)
      if(o2cal@var != "") SumTable = fitSumTable(o2cal, "O2", SumTable)
      n2cal = getCal("N2", SumTable)
      if(n2cal@var != "") SumTable = fitSumTable(n2cal, "N2", SumTable)
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Now write SumTable to .csv file
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      outfile = choose.files(mydir@dO, filters = myFilters, caption = "Choose SPKYNetFit output file")
      write.csv(SumTable, file = outfile, quote = FALSE, row.names = FALSE)
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Create reduced summary table with only sample info and amounts
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      tiny.filename = strsplit(outfile, ".csv")
      tiny.filename = paste(tiny.filename, "_tiny.csv", sep = "") 
      if (gc == 3){
        column.list = c("File Name", "DateTime", "Sample Name", "Sample Number")
        tiny.SumTable = cbind(SumTable[,column.list],SumTable[,grep("Amount",colnames(SumTable), ignore.case = FALSE)])    
      }else{
        column.list = c("Run", "File Name", "Inj. Date/Time", "Sample Name")
        tiny.SumTable = cbind(SumTable[,column.list],SumTable[,grep("Amount",colnames(SumTable), ignore.case = FALSE)])
      }
      write.csv(tiny.SumTable, file = tiny.filename, quote = FALSE, row.names = FALSE)
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Now extract and write STD data to SparkySTD.csv file
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      stdfn = file(paste(mydir@dP, "/SparkyJrSTD.csv", sep = ""), "a")
      write.csv(SumTable[which(stds),], file = stdfn, quote = FALSE, row.names = FALSE)
      close(stdfn)
      cat("\nSparky Jr. parsed data written to: \n\n", outfile, "\n\n")
      cat("\nParsed by ", computer, "\n")
      
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Save parsing info to parserlog.txt and parsed.txt in ECDData folder
      # Save parsing info to log file. One file for all parsers?
      # Info:  date, parser name, sequence name, output path/file name, errors
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      fn = file(paste(mydir@dP, "/parserlog.txt", sep = ""), "a")
      sequenceName = strsplit(infile2, "\\\\")[[1]][length(strsplit(infile2, "\\\\")[[1]]) - 1]
      writeLines(paste(date(), systemName, sequenceName, outfile, parserError, sep = ","), fn)
      close(fn)
      bname = substr(basename(infile2), 1, nchar(basename(infile2)))
      fn = file(paste(substr(infile2,1, nchar(infile2) - nchar(bname)), "parsed.txt", sep = ""))
      writeLines(paste(parserVersion, date(), systemName, sequenceName, outfile, computer, parserError, sep = ","), fn)
      close(fn)
      
      # End of Recal portion of script.
      
      
    }else{
      tkmessageBox(title="Parser error!!", message = "No data were processed!  Fix problem with data files then try parser again.", icon="warning", type="ok")
      tkmessageBox(title="HPCHEM-only output.", message = "You can output HPCHEM data, but no calculations will be performed.", icon="warning", type="ok")
      
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      # Now write partial SumTable to .csv file (HPChem data only)
      #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      outfile = choose.files(mydir@dO, filters = myFilters, caption = "Choose HPCHEM-only output file")
      write.csv(SumTable, file = outfile, quote = FALSE, row.names = FALSE)
    }
    
  }else{
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    # Manual calibration selected. Output raw data only, to be processed by hand.
    # Now write SumTable to .csv file
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    
    HTML.title("Data required manual calibration.", file = HTMLoutput)
    
    outfile = choose.files(mydir@dO, filters = myFilters, caption = "Choose SPKYNetFit output file")
    write.csv(SumTable, file = outfile, quote = FALSE, row.names = FALSE)
    
  }

#End of data processing

setwd(mydir@dW)