# Lachat Drift Correction Codebook

The R script reads a user selected raw data file (.csv) generated by Onmion version 3.0.  The raw data file will contain data for 1 to 3 of the following compounds:  NH3, NO3 or PO4. 

####Table 1.
|  column index |||Variable names|||
|---|---|---|---|---|---|
|[1]| Sample.ID| Replicate.Number |Cup.Number  |Detection.Date| Detection.Time  |     
|[6]|Auto.Dilution.Factor|Analyte.Name|Peak.Area|   Peak.Height|Peak.Concentration|   
|[11]|Concentration.Units|  Channel.Number| 
**Table 1. Lists the variable names in the raw data file.**


Analyte.Name, Peak.Area, Peak.Height, Peak.Concentration, Concentration.Units, and Channel.Number are repeated for each compound and will have a *.1* or *.2* appended for each additional analyte.

The script renames the columns to include the analyte name and removes the "." from all the names.  It then reformats the data into tidy form so that it only contains 6 analyte related variables:

####Table 2.
|  column index |||Variable names|||
|---|---|---|---|---|---|
|[1]| SampleID| ReplicateNumber |CupNumber  |DetectionDate| DetectionTime  |     
|[6]|AutoDilutionFactor|Area|Height|Concentration|Units|
|[11]|Channel|Analyte| 
**Table 2. Lists the variable names in the tidy data file.**

The Lachat can sometimes exhibit response drift over time.  If the standard error for the concentration of the highest standard is strictly larger than 10%, the concentrations for the analyte are corrected using a polynomial fit scale factor. **This correction is still under development.**

The script exports 3 data files.  Each file has the name of the original data file plus a suffix to identify the content of the file:

_tidy :  A tidy data file.
  
_messy :  A messy data file that saves the data in the original layout but with renamed columns discussed above.
 
_LDC :  A drift corrected data file.


---   
Since the drift correction is still under development the _LDC file is a copy of the tidy data file.

