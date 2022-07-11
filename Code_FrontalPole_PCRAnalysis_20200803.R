#Adriana's Frontal Pole PCR Data
## Goal: Figure out why the "housekeeping genes" are showing large amounts of variability
## Megan Hagenauer
## August 3, 2020


####################

#Reading in Files:

##This was re-run 2020 08 25 because I discovered that we had been missing data from card 22.eds (sample IDs 41-44) on the server. Adriana reoutputted it from the proprietary website and uploaded it, so it is now included.

#Reading in general experiment info:

setwd("~/Documents/Microarray Gen/FrontalPole")

HousekeepingGenes<-read.csv("HouseKeepingGenes_UsedInPCR_FrontalPole.csv", header=T, stringsAsFactors = F)
str(HousekeepingGenes)

# 'data.frame':	24 obs. of  1 variable:
#   $ GeneSymbol: chr  "18S" "ACTB" "B2M" "CASC3" ...

SubjectInfo<-read.csv("Balancing RNA Extraction04162018_UseThis_Simplified.csv", header=T, stringsAsFactors = F)
str(SubjectInfo)

# 'data.frame':	72 obs. of  12 variables:
#   $ Barcode        : chr  "B007375A" "B010169A" "B000014A" "B010842A" ...
# $ Subject.Number : int  5000 5313 2169 4087 5066 2950 4383 4619 4819 3711 ...
# $ Cohort         : chr  "Cohort 11" "Cohort 13" "Dep Cohort 1" "Dep Cohort 6" ...
# $ Diagnosis      : chr  "Control" "Schiz" "Control" "BP" ...
# $ Age            : int  77 42 18 51 47 45 40 34 79 59 ...
# $ Gender         : chr  "M" "M" "M" "F" ...
# $ pH             : num  6.79 6.77 6.97 6.81 6.6 7.05 6.77 6.7 6.75 6.55 ...
# $ AFS            : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final    : num  16.5 23.7 22 16.8 21 20 26.3 17.2 21 27 ...
# $ Slab.Format    : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number    : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group: int  1 1 1 1 1 1 2 2 2 2 ...

SubjectInfo$Diagnosis<-as.factor(SubjectInfo$Diagnosis)
levels(SubjectInfo$Diagnosis)
SubjectInfo$Diagnosis<-relevel(SubjectInfo$Diagnosis, ref="Control")
levels(SubjectInfo$Diagnosis)
#[1] "Control" "BP"      "Schiz" 


#Decoders were made with the following info:

#The original decoding file ("Barcode Order Form.xlsx") needed to be reformatted:
#From Adriana: Cards 1-36 are the first set (GABA glutamate) and have 4 samples per card, then the duplicate for each sample is in the next card. For example sample 1 is on card 1 and card 2.
#Cards 39 to 58 are the Dopamine/Serotonin, those have 8 samples per card then the duplicates for those samples are in the following card. There are no cards 37 or 38. 
#Also, samples 11, 35, and 36 did not work on the extraction RNA process so they were redone and relabeled as samples 73, 74, and 75.

#In the folder for each card, the raw CT values are in the spreadsheet called Results. The wells that did not amplify would say No Amp or Undetermined instead of a CT number. 
#Samples 27 and 38 were the ones that had a whole stem of wells not amplified (meaning one of the duplicates didn't work). 


Decoder_GabaGluCards<-read.csv("Decoder_FrontalPole_GabaGlu_Cards.csv", header=T, stringsAsFactors = F)
str(Decoder_GabaGluCards)
#' #'data.frame':	144 obs. of  5 variables:
#'  $ Subject.Number       : int  5000 5000 5313 5313 2169 2169 4087 4087 5066 5066 ...
#' $ ID                   : int  1 1 2 2 3 3 4 4 5 5 ...
#' $ Folder               : chr  "Card 1.eds" "Card 2.eds" "Card 1.eds" "Card 2.eds" ...
#' $ SampleNumber         : int  1 1 2 2 3 3 4 4 1 1 ...
#' $ QC_AmplificationIssue: chr  "N" "N" "N" "N" ..."

Decoder_DA5HTCards<-read.csv("Decoder_FrontalPole_DA5HT_Cards.csv", header=T, stringsAsFactors = F)
str(Decoder_DA5HTCards)
# 'data.frame':	144 obs. of  5 variables:
#   $ Subject.Number       : int  5000 5000 5313 5313 2169 2169 4087 4087 5066 5066 ...
# $ ID                   : int  1 1 2 2 3 3 4 4 5 5 ...
# $ Folder               : logi  NA NA NA NA NA NA ...
# $ SampleNumber         : logi  NA NA NA NA NA NA ...
# $ QC_AmplificationIssue: chr  "N" "N" "N" "N" ...

#Note: I don't have the folder entered for this one or sample number because it didn't follow the same predictable order as the GabaGlu Cards (which is probably good)


#***********************************


#Reading in the individual results files for the GabaGlu dataset:


#We need to loop over all folders in the directory and grab the results file from within each folder.

setwd("~/Documents/Microarray Gen/FrontalPole/PCR")

foldernames<-cbind(list.files(), list.files(full.names=T))
str(foldernames)
#chr [1:56, 1:2] "Card 1.eds" "Card 10.eds" "Card 11.eds" "Card 12.eds" "Card 13.eds" "Card 14.eds" "Card 15.eds" ...
head(foldernames)
#          [,1]          [,2]           
# [1,] "Card 1.eds"  "./Card 1.eds" 
# [2,] "Card 10.eds" "./Card 10.eds"
# [3,] "Card 11.eds" "./Card 11.eds"
# [4,] "Card 12.eds" "./Card 12.eds"
# [5,] "Card 13.eds" "./Card 13.eds"
# [6,] "Card 14.eds" "./Card 14.eds"

setwd(foldernames[1,2])
list.files()
#[1] "Amplification Data.csv" "Plate QC.csv"           "Results.csv"  

setwd("~/Documents/Microarray Gen/FrontalPole")
write.csv(foldernames, "PCR_foldernames.csv")
#I quickly annotated this with the biological pathway targeted by each card, and then read it back in:

foldernames<-read.csv("PCR_foldernames.csv", header=T, stringsAsFactors = F)
str(foldernames)

# 'data.frame':	56 obs. of  4 variables:
#   $ X              : int  1 2 3 4 5 6 7 8 9 10 ...
# $ FileName       : chr  "Card 1.eds" "Card 10.eds" "Card 11.eds" "Card 12.eds" ...
# $ FolderPath     : chr  "./Card 1.eds" "./Card 10.eds" "./Card 11.eds" "./Card 12.eds" ...
# $ TargetedBiology: chr  "GabaGlu" "GabaGlu" "GabaGlu" "GabaGlu" ...

foldernames_GabaGlu<-foldernames[foldernames$TargetedBiology=="GabaGlu",]
foldernames_DA5HT<-foldernames[foldernames$TargetedBiology=="DA5HT",]

setwd("~/Documents/Microarray Gen/FrontalPole/PCR")
setwd(foldernames_GabaGlu$FolderPath[1])

#We need to ignore the first 5 rows when reading in data.
temp<-read.table("Results.csv", sep=",", skip=5, header=T, stringsAsFactors = F)
str(temp)

# 'data.frame':	384 obs. of  23 variables:
#   $ Well                       : chr  "A1" "A2" "A3" "A4" ...
# $ Sample.Name                : chr  "No Sample" "1" "1" "1" ...
# $ Target.Name                : chr  "ABAT-Hs00609436_m1" "ADCY7-Hs00936808_m1" "ADORA1-Hs00379752_m1" "ADORA2A-Hs00169123_m1" ...
# $ Amp.Score                  : num  1.28 1.25 1.14 1.23 1.23 ...
# $ Amp.Status                 : chr  "Amp" "Amp" "Amp" "Amp" ...
# $ Task                       : chr  "UNKNOWN" "UNKNOWN" "UNKNOWN" "UNKNOWN" ...
# $ Cq                         : chr  "20.513" "25.734" "23.737" "27.270" ...
# $ Cq.Mean                    : chr  "20.513" "25.734" "23.737" "27.27" ...
# $ Cq.Standard.Deviation      : chr  "-" "-" "-" "-" ...
# $ Quantity                   : chr  "-" "-" "-" "-" ...
# $ Quantity.Mean              : chr  "-" "-" "-" "-" ...
# $ Quantity.Standard.Deviation: chr  "-" "-" "-" "-" ...
# $ Intercept                  : chr  "-" "-" "-" "-" ...
# $ R.Squared                  : chr  "-" "-" "-" "-" ...
# $ Slope                      : chr  "-" "-" "-" "-" ...
# $ Efficiency                 : chr  "-" "-" "-" "-" ...
# $ Auto.Threshold             : chr  "true" "true" "true" "true" ...
# $ Threshold                  : num  0.183 0.161 0.09 0.175 0.133 0.193 0.14 0.106 0.195 0.136 ...
# $ Auto.Baseline              : chr  "true" "true" "true" "true" ...
# $ Baseline.Start             : num  3 3 3 3 3 3 3 3 3 3 ...
# $ Baseline.End               : num  15 20 19 22 18 17 19 19 16 19 ...
# $ Omit                       : chr  "false" "false" "false" "false" ...
# $ X                          : logi  NA NA NA NA NA NA ...

temp$Card<-rep(foldernames_GabaGlu$FileName[1], length(temp[,1]))
 
str(temp) 

#Let's use that as the start for the df:
Concatenated_GabaGlu<-temp
rm(temp)  
#Cool, that seems to work.


for(i in c(2:length(foldernames_GabaGlu$FileName))){
  setwd("~/Documents/Microarray Gen/FrontalPole/PCR")
  setwd(foldernames_GabaGlu$FolderPath[i])
  print(getwd())
  
  temp<-read.table("Results.csv", sep=",", skip=5, header=T, stringsAsFactors = F)
  temp$Card<-rep(foldernames_GabaGlu$FileName[i], length(temp[,1]))
  
  Concatenated_GabaGlu<-rbind.data.frame(Concatenated_GabaGlu, temp)
  rm(temp)
  
}

# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 10.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 11.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 12.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 13.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 14.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 15.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 16.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 17.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 18.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 19.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 2.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 20.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 21.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 22.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 23.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 24.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 25.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 26.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 27.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 28.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 29.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 3.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 30.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 31.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 32.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 33.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 34.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 35.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 36.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 4.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 5.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 6.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 7.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 8.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 9.eds"

str(Concatenated_GabaGlu)

#'data.frame':	13824 obs. of  24 variables:
#   $ Well                       : chr  "A1" "A2" "A3" "A4" ...
# $ Sample.Name                : chr  "No Sample" "1" "1" "1" ...
# $ Target.Name                : chr  "ABAT-Hs00609436_m1" "ADCY7-Hs00936808_m1" "ADORA1-Hs00379752_m1" "ADORA2A-Hs00169123_m1" ...
# $ Amp.Score                  : num  1.28 1.25 1.14 1.23 1.23 ...
# $ Amp.Status                 : chr  "Amp" "Amp" "Amp" "Amp" ...
# $ Task                       : chr  "UNKNOWN" "UNKNOWN" "UNKNOWN" "UNKNOWN" ...
# $ Cq                         : chr  "20.513" "25.734" "23.737" "27.270" ...
# $ Cq.Mean                    : chr  "20.513" "25.734" "23.737" "27.27" ...
# $ Cq.Standard.Deviation      : chr  "-" "-" "-" "-" ...
# $ Quantity                   : chr  "-" "-" "-" "-" ...
# $ Quantity.Mean              : chr  "-" "-" "-" "-" ...
# $ Quantity.Standard.Deviation: chr  "-" "-" "-" "-" ...
# $ Intercept                  : chr  "-" "-" "-" "-" ...
# $ R.Squared                  : chr  "-" "-" "-" "-" ...
# $ Slope                      : chr  "-" "-" "-" "-" ...
# $ Efficiency                 : chr  "-" "-" "-" "-" ...
# $ Auto.Threshold             : chr  "true" "true" "true" "true" ...
# $ Threshold                  : num  0.183 0.161 0.09 0.175 0.133 0.193 0.14 0.106 0.195 0.136 ...
# $ Auto.Baseline              : chr  "true" "true" "true" "true" ...
# $ Baseline.Start             : num  3 3 3 3 3 3 3 3 3 3 ...
# $ Baseline.End               : num  15 20 19 22 18 17 19 19 16 19 ...
# $ Omit                       : chr  "false" "false" "false" "false" ...
# $ X                          : logi  NA NA NA NA NA NA ...
# $ Card                       : chr  "Card 1.eds" "Card 1.eds" "Card 1.eds" "Card 1.eds" ...


#Doing some quick sanity checks to make sure everything read in properly:

table(Concatenated_GabaGlu$Card)

# Card 1.eds Card 10.eds Card 11.eds Card 12.eds Card 13.eds Card 14.eds Card 15.eds Card 16.eds Card 17.eds Card 18.eds 
# 384         384         384         384         384         384         384         384         384         384 
# Card 19.eds  Card 2.eds Card 20.eds Card 21.eds Card 22.eds Card 23.eds Card 24.eds Card 25.eds Card 26.eds Card 27.eds 
# 384         384         384         384         384         384         384         384         384         384 
# Card 28.eds Card 29.eds  Card 3.eds Card 30.eds Card 31.eds Card 32.eds Card 33.eds Card 34.eds Card 35.eds Card 36.eds 
# 384         384         384         384         384         384         384         384         384         384 
# Card 4.eds  Card 5.eds  Card 6.eds  Card 7.eds  Card 8.eds  Card 9.eds 
# 384         384         384         384         384         384 

table(Concatenated_GabaGlu$Sample.Name)

# 1        10        12        13        14        15        16        17        18        19         2        20 
# 190       192       192       190       192       192       192       190       192       192       192       192 
# 21        22        23        24        25        26        27        28        29         3        30        31 
# 190       192       192       192       190       192       192       192       190       192       192       192 
# 32        33        34        37        38        39         4        40        41        42        43        44 
# 192       190       192       190       192       192       192       192       190       192       192       192 
# 45        46        47        48        49         5        50        51        52        53        54        55 
# 190       192       192       192       190       190       192       192       192       190       192       192 
# 56        57        58        59         6        60        61        62        63        64        65        66 
# 192       190       192       192       192       192       190       192       192       192       190       192 
# 67        68        69         7        70        71        72        73        74        75         8         9 
# 192       192       190       192       192       192       192       192       192       192       192       190 
# No Sample 
# 36 

#Hmmm... I wonder why there are different #'s of rows of data for different samples (190-192 samples)


HowManyMeasurementsPerSample<-table(Concatenated_GabaGlu$Card, Concatenated_GabaGlu$Sample.Name)

HowManyMeasurementsPerSample

# 1 10 12 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31 32 33 34 37 38 39  4 40 41 42 43 44
# Card 1.eds  95  0  0  0  0  0  0  0  0  0 96  0  0  0  0  0  0  0  0  0  0 96  0  0  0  0  0  0  0  0 96  0  0  0  0  0
# Card 10.eds  0  0  0  0  0  0  0 95 96 96  0 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 11.eds  0  0  0  0  0  0  0  0  0  0  0  0 95 96 96 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 12.eds  0  0  0  0  0  0  0  0  0  0  0  0 95 96 96 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 13.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 95 96 96 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 14.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 95 96 96 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 15.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 95  0 96 96 96  0  0  0  0  0  0  0  0  0  0  0
# Card 16.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 95  0 96 96 96  0  0  0  0  0  0  0  0  0  0  0
# Card 17.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 95 96  0  0  0  0  0  0  0  0  0
# Card 18.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 95 96  0  0  0  0  0  0  0  0  0
# Card 19.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 95 96 96  0 96  0  0  0  0
# Card 2.eds  95  0  0  0  0  0  0  0  0  0 96  0  0  0  0  0  0  0  0  0  0 96  0  0  0  0  0  0  0  0 96  0  0  0  0  0
# Card 20.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 95 96 96  0 96  0  0  0  0
# 
# 45 46 47 48 49  5 50 51 52 53 54 55 56 57 58 59  6 60 61 62 63 64 65 66 67 68 69  7 70 71 72 73 74 75  8  9
# Card 1.eds   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 10.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 11.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 12.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 13.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 14.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 15.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 16.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 17.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 96 96  0  0
# Card 18.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 96 96  0  0
# Card 19.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 2.eds   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 20.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# 
# No Sample
# Card 1.eds          1
# Card 10.eds         1
# Card 11.eds         1
# Card 12.eds         1
# Card 13.eds         1
# Card 14.eds         1
# Card 15.eds         1
# Card 16.eds         1
# Card 17.eds         1
# Card 18.eds         1
# Card 19.eds         1
# Card 2.eds          1
# Card 20.eds         1
# [ reached getOption("max.print") -- omitted 23 rows ]

setwd("~/Documents/Microarray Gen/FrontalPole")
write.csv(HowManyMeasurementsPerSample, "HowManyMeasurementsPerSample.csv")

#Looks like some samples have 95 measurements and some have 96. Huh. 
#perhaps the header is being included in the rbind?
#Alright, it looks like the first sample on each card is missing 1 row of data. Maybe I have skip set to high.
#Nope - the problem is in the original results spreadsheets. The first row of data is always labeled "No Sample" but appears to actually be the measurement for gene "ABAT" for the first sample in the dataset.
#I asked Adriana about it - she says that it is a labeling error... so how do we fix it easily?  Let's grab the label from one row beneath it...

sum(Concatenated_GabaGlu$Sample.Name=="No Sample")
#[1] 36

which(Concatenated_GabaGlu$Sample.Name=="No Sample")
which(Concatenated_GabaGlu$Sample.Name=="No Sample")+1

IndicesForCorrectedSampleIDs<-which(Concatenated_GabaGlu$Sample.Name=="No Sample")+1

Concatenated_GabaGlu$Sample.Name[IndicesForCorrectedSampleIDs]

Concatenated_GabaGlu$Sample.Name[Concatenated_GabaGlu$Sample.Name=="No Sample"]

Concatenated_GabaGlu$Sample.Name[Concatenated_GabaGlu$Sample.Name=="No Sample"]<-Concatenated_GabaGlu$Sample.Name[IndicesForCorrectedSampleIDs]
  
table(Concatenated_GabaGlu$Sample.Name)

# 1  10  12  13  14  15  16  17  18  19   2  20  21  22  23  24  25  26  27  28  29   3  30  31  32  33  34  37  38  39 
# 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 
# 4  40  41  42  43  44  45  46  47  48  49   5  50  51  52  53  54  55  56  57  58  59   6  60  61  62  63  64  65  66 
# 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 192 
# 67  68  69   7  70  71  72  73  74  75   8   9 
# 192 192 192 192 192 192 192 192 192 192 192 192 

#Fixed. :)

HowManyMeasurementsPerSample<-table(Concatenated_GabaGlu$Card, Concatenated_GabaGlu$Sample.Name)

setwd("~/Documents/Microarray Gen/FrontalPole")
write.csv(HowManyMeasurementsPerSample, "HowManyMeasurementsPerSample_GabaGlu_Fixed.csv")


#############

#Reading in the individual files for the DA5HT dataset:


setwd("~/Documents/Microarray Gen/FrontalPole/PCR")
setwd(foldernames_DA5HT$FolderPath[1])

#We need to ignore the first 5 rows when reading in data.
temp<-read.table("Results.csv", sep=",", skip=5, header=T, stringsAsFactors = F)
str(temp)


# 'data.frame':	384 obs. of  23 variables:
#   $ Well                       : chr  "A1" "A2" "A3" "A4" ...
# $ Sample.Name                : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Target.Name                : chr  "ADRB1-Hs02330048_s1" "ADRB2-Hs00240532_s1" "COMT-Hs00241349_m1" "DBH-Hs01089840_m1" ...
# $ Amp.Score                  : num  1.25 1.38 1.09 1.31 1.42 ...
# $ Amp.Status                 : chr  "Amp" "Amp" "Amp" "Amp" ...
# $ Task                       : chr  "UNKNOWN" "UNKNOWN" "UNKNOWN" "UNKNOWN" ...
# $ Cq                         : chr  "24.245" "25.087" "24.579" "27.823" ...
# $ Cq.Mean                    : chr  "24.245" "25.087" "24.579" "27.823" ...
# $ Cq.Standard.Deviation      : chr  "-" "-" "-" "-" ...
# $ Quantity                   : chr  "-" "-" "-" "-" ...
# $ Quantity.Mean              : chr  "-" "-" "-" "-" ...
# $ Quantity.Standard.Deviation: chr  "-" "-" "-" "-" ...
# $ Intercept                  : chr  "-" "-" "-" "-" ...
# $ R.Squared                  : chr  "-" "-" "-" "-" ...
# $ Slope                      : chr  "-" "-" "-" "-" ...
# $ Efficiency                 : chr  "-" "-" "-" "-" ...
# $ Auto.Threshold             : chr  "true" "true" "true" "true" ...
# $ Threshold                  : num  0.163 0.239 0.048 0.167 0.24 0.202 0.227 0.194 0.184 0.228 ...
# $ Auto.Baseline              : chr  "true" "true" "true" "true" ...
# $ Baseline.Start             : num  3 3 3 3 3 3 3 3 3 3 ...
# $ Baseline.End               : num  19 20 21 23 27 20 23 27 27 21 ...
# $ Omit                       : chr  "false" "false" "false" "false" ...
# $ X                          : logi  NA NA NA NA NA NA ...

temp$Card<-rep(foldernames_DA5HT$FileName[1], length(temp[,1]))

str(temp) 

#Let's use that as the start for the df:
Concatenated_DA5HT<-temp
rm(temp)  
#Cool, that seems to work.


for(i in c(2:length(foldernames_DA5HT$FileName))){
  setwd("~/Documents/Microarray Gen/FrontalPole/PCR")
  setwd(foldernames_DA5HT$FolderPath[i])
  print(getwd())
  
  temp<-read.table("Results.csv", sep=",", skip=5, header=T, stringsAsFactors = F)
  temp$Card<-rep(foldernames_DA5HT$FileName[i], length(temp[,1]))
  
  Concatenated_DA5HT<-rbind.data.frame(Concatenated_DA5HT, temp)
  rm(temp)
  
}

# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 40.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 41.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 42.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 43.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 44.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 45.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 46.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 47.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 48.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 49.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 50.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 51.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 52.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 53.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 54.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 55.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 56.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 57.eds"
# [1] "/Users/mhh/Documents/Microarray Gen/FrontalPole/PCR/Card 58.eds"


str( Concatenated_DA5HT)

# 'data.frame':	7680 obs. of  24 variables:
#   $ Well                       : chr  "A1" "A2" "A3" "A4" ...
# $ Sample.Name                : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Target.Name                : chr  "ADRB1-Hs02330048_s1" "ADRB2-Hs00240532_s1" "COMT-Hs00241349_m1" "DBH-Hs01089840_m1" ...
# $ Amp.Score                  : num  1.25 1.38 1.09 1.31 1.42 ...
# $ Amp.Status                 : chr  "Amp" "Amp" "Amp" "Amp" ...
# $ Task                       : chr  "UNKNOWN" "UNKNOWN" "UNKNOWN" "UNKNOWN" ...
# $ Cq                         : chr  "24.245" "25.087" "24.579" "27.823" ...
# $ Cq.Mean                    : chr  "24.245" "25.087" "24.579" "27.823" ...
# $ Cq.Standard.Deviation      : chr  "-" "-" "-" "-" ...
# $ Quantity                   : chr  "-" "-" "-" "-" ...
# $ Quantity.Mean              : chr  "-" "-" "-" "-" ...
# $ Quantity.Standard.Deviation: chr  "-" "-" "-" "-" ...
# $ Intercept                  : chr  "-" "-" "-" "-" ...
# $ R.Squared                  : chr  "-" "-" "-" "-" ...
# $ Slope                      : chr  "-" "-" "-" "-" ...
# $ Efficiency                 : chr  "-" "-" "-" "-" ...
# $ Auto.Threshold             : chr  "true" "true" "true" "true" ...
# $ Threshold                  : num  0.163 0.239 0.048 0.167 0.24 0.202 0.227 0.194 0.184 0.228 ...
# $ Auto.Baseline              : chr  "true" "true" "true" "true" ...
# $ Baseline.Start             : num  3 3 3 3 3 3 3 3 3 3 ...
# $ Baseline.End               : num  19 20 21 23 27 20 23 27 27 21 ...
# $ Omit                       : chr  "false" "false" "false" "false" ...
# $ X                          : logi  NA NA NA NA NA NA ...
# $ Card                       : chr  "Card 39.eds" "Card 39.eds" "Card 39.eds" "Card 39.eds" ...


#Doing some quick sanity checks:

table(Concatenated_DA5HT$Card)

# Card 39.eds Card 40.eds Card 41.eds Card 42.eds Card 43.eds Card 44.eds Card 45.eds Card 46.eds Card 47.eds Card 48.eds 
# 384         384         384         384         384         384         384         384         384         384 
# Card 49.eds Card 50.eds Card 51.eds Card 52.eds Card 53.eds Card 54.eds Card 55.eds Card 56.eds Card 57.eds Card 58.eds 
# 384         384         384         384         384         384         384         384         384         384 


table(Concatenated_DA5HT$Sample.Name)

# 1   2   3   4   5   6   7   8   9  10  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33 
# 192  96  96  96  96  96  96  96  96  96  96  96  96  96  96  96 192  96  96  96 192  96  96  96  96  96  96  96  96 192  96  96 
# 34  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67 
# 96  96  96  96  96 192  96  96  96  96 192  96  96  96  96  96  96  96  96  96  96  96 192  96  96 192  96  96  96  96  96  96 
# 68  69  70  71  72  73  74  75 
# 96  96  96  96  96  96  96  96 

#Several samples have 192 measurements: sample numbers 1, 18, 22, 31, 41, 46, 58, 61
#Why...?

table(Concatenated_DA5HT$Card, Concatenated_DA5HT$Sample.Name)

#             1  2  3  4  5  6  7  8  9 10 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 37 38 39 40 41
# Card 39.eds 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 40.eds 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 41.eds  0  0  0  0  0  0  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 42.eds  0  0  0  0  0  0  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 43.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 44.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 45.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0
# Card 46.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0
# Card 47.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 48 48 48 48 48 48
# Card 48.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 48 48 48 48 48 48
# Card 49.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 50.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 51.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# 
#             42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
# Card 39.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 40.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 41.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 42.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 43.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 44.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 45.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 46.eds  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 47.eds 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 48.eds 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 49.eds  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 50.eds  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# Card 51.eds  0  0  0  0  0  0  0  0  0  0 48 48 48 48 48 48 48 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

#...


HowManyMeasurementsPerSample_DA5HT<-table(Concatenated_DA5HT$Card, Concatenated_DA5HT$Sample.Name)

setwd("~/Documents/Microarray Gen/FrontalPole")
write.csv(HowManyMeasurementsPerSample_DA5HT, "HowManyMeasurementsPerSample_DA5HT.csv")

#most samples were run in duplicate (2 cards), but 8 samples (sample numbers 1, 18, 22, 31, 41, 46, 58, 61) were run in quadruplicate (4 cards), with cards 57 and 58 entirely composed of samples that had been run on other cards.

#I asked Adriana, and she says that she had two extra cards so she re-ran any samples that seemed like there could have been issues with pippetting/well-filling (without having seen the results). So hypothetically we could either use the later cards for those samples, or use all 4 quadruplicates if the data looks clean.


############################

#We need to take the name of the sample number and join it with the subject info:


#First, I changed the Sample.Name column to "ID" to make it match the sample info data.frame for later joining:
str(Concatenated_GabaGlu)
colnames(Concatenated_GabaGlu)[2]<-"ID"

str(Concatenated_DA5HT)
colnames(Concatenated_DA5HT)[2]<-"ID"


library(plyr)

str(SubjectInfo)

# 'data.frame':	72 obs. of  12 variables:
#   $ Barcode        : chr  "B007375A" "B010169A" "B000014A" "B010842A" ...
# $ Subject.Number : int  5000 5313 2169 4087 5066 2950 4383 4619 4819 3711 ...
# $ Cohort         : chr  "Cohort 11" "Cohort 13" "Dep Cohort 1" "Dep Cohort 6" ...
# $ Diagnosis      : chr  "Control" "Schiz" "Control" "BP" ...
# $ Age            : int  77 42 18 51 47 45 40 34 79 59 ...
# $ Gender         : chr  "M" "M" "M" "F" ...
# $ pH             : num  6.79 6.77 6.97 6.81 6.6 7.05 6.77 6.7 6.75 6.55 ...
# $ AFS            : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final    : num  16.5 23.7 22 16.8 21 20 26.3 17.2 21 27 ...
# $ Slab.Format    : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number    : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group: int  1 1 1 1 1 1 2 2 2 2 ...

str(Decoder_GabaGluCards)

# 'data.frame':	144 obs. of  5 variables:
#   $ Subject.Number       : int  5000 5000 5313 5313 2169 2169 4087 4087 5066 5066 ...
# $ ID                   : int  1 1 2 2 3 3 4 4 5 5 ...
# $ Folder               : chr  "Card 1.eds" "Card 2.eds" "Card 1.eds" "Card 2.eds" ...
# $ SampleNumber         : int  1 1 2 2 3 3 4 4 1 1 ...
# $ QC_AmplificationIssue: chr  "N" "N" "N" "N" ...

#Twice the size because there are two rows for each subject.


Decoder_GabaGluCards_w_SubjectInfo<-join(SubjectInfo, Decoder_GabaGluCards[,c(1:2,4:5)], by="Subject.Number", type="left", match="first")
str(Decoder_GabaGluCards_w_SubjectInfo)

# 'data.frame':	72 obs. of  15 variables:
#   $ Subject.Number       : int  5000 5313 2169 4087 5066 2950 4383 4619 4819 3711 ...
# $ Barcode              : chr  "B007375A" "B010169A" "B000014A" "B010842A" ...
# $ Cohort               : chr  "Cohort 11" "Cohort 13" "Dep Cohort 1" "Dep Cohort 6" ...
# $ Diagnosis            : chr  "Control" "Schiz" "Control" "BP" ...
# $ Age                  : int  77 42 18 51 47 45 40 34 79 59 ...
# $ Gender               : chr  "M" "M" "M" "F" ...
# $ pH                   : num  6.79 6.77 6.97 6.81 6.6 7.05 6.77 6.7 6.75 6.55 ...
# $ AFS                  : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final          : num  16.5 23.7 22 16.8 21 20 26.3 17.2 21 27 ...
# $ Slab.Format          : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number          : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group      : int  1 1 1 1 1 1 2 2 2 2 ...
# $ ID                   : int  1 2 3 4 5 6 7 8 9 10 ...
# $ SampleNumber         : int  1 2 3 4 1 2 3 4 1 2 ...
# $ QC_AmplificationIssue: chr  "N" "N" "N" "N" ...


write.csv(Decoder_GabaGluCards_w_SubjectInfo, "Decoder_GabaGluCards_w_SubjectInfo.csv")


Concatenated_GabaGlu_wSubjectInfo<-join(Concatenated_GabaGlu, Decoder_GabaGluCards_w_SubjectInfo, by="ID", type="left")
str(Concatenated_GabaGlu_wSubjectInfo)

# 'data.frame':	13824 obs. of  38 variables:
#   $ Well                       : chr  "A1" "A2" "A3" "A4" ...
# $ ID                         : chr  "1" "1" "1" "1" ...
# $ Target.Name                : chr  "ABAT-Hs00609436_m1" "ADCY7-Hs00936808_m1" "ADORA1-Hs00379752_m1" "ADORA2A-Hs00169123_m1" ...
# $ Amp.Score                  : num  1.28 1.25 1.14 1.23 1.23 ...
# $ Amp.Status                 : chr  "Amp" "Amp" "Amp" "Amp" ...
# $ Task                       : chr  "UNKNOWN" "UNKNOWN" "UNKNOWN" "UNKNOWN" ...
# $ Cq                         : chr  "20.513" "25.734" "23.737" "27.270" ...
# $ Cq.Mean                    : chr  "20.513" "25.734" "23.737" "27.27" ...
# $ Cq.Standard.Deviation      : chr  "-" "-" "-" "-" ...
# $ Quantity                   : chr  "-" "-" "-" "-" ...
# $ Quantity.Mean              : chr  "-" "-" "-" "-" ...
# $ Quantity.Standard.Deviation: chr  "-" "-" "-" "-" ...
# $ Intercept                  : chr  "-" "-" "-" "-" ...
# $ R.Squared                  : chr  "-" "-" "-" "-" ...
# $ Slope                      : chr  "-" "-" "-" "-" ...
# $ Efficiency                 : chr  "-" "-" "-" "-" ...
# $ Auto.Threshold             : chr  "true" "true" "true" "true" ...
# $ Threshold                  : num  0.183 0.161 0.09 0.175 0.133 0.193 0.14 0.106 0.195 0.136 ...
# $ Auto.Baseline              : chr  "true" "true" "true" "true" ...
# $ Baseline.Start             : num  3 3 3 3 3 3 3 3 3 3 ...
# $ Baseline.End               : num  15 20 19 22 18 17 19 19 16 19 ...
# $ Omit                       : chr  "false" "false" "false" "false" ...
# $ X                          : logi  NA NA NA NA NA NA ...
# $ Card                       : chr  "Card 1.eds" "Card 1.eds" "Card 1.eds" "Card 1.eds" ...
# $ Subject.Number             : int  5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 ...
# $ Barcode                    : chr  "B007375A" "B007375A" "B007375A" "B007375A" ...
# $ Cohort                     : chr  "Cohort 11" "Cohort 11" "Cohort 11" "Cohort 11" ...
# $ Diagnosis                  : chr  "Control" "Control" "Control" "Control" ...
# $ Age                        : int  77 77 77 77 77 77 77 77 77 77 ...
# $ Gender                     : chr  "M" "M" "M" "M" ...
# $ pH                         : num  6.79 6.79 6.79 6.79 6.79 6.79 6.79 6.79 6.79 6.79 ...
# $ AFS                        : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final                : num  16.5 16.5 16.5 16.5 16.5 16.5 16.5 16.5 16.5 16.5 ...
# $ Slab.Format                : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number                : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group            : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SampleNumber               : int  1 1 1 1 1 1 1 1 1 1 ...
# $ QC_AmplificationIssue      : chr  "N" "N" "N" "N" ...



Decoder_DA5HTCards_w_SubjectInfo<-join(SubjectInfo, Decoder_DA5HTCards[,c(1:2,4:5)], by="Subject.Number", type="left", match="first")
str(Decoder_DA5HTCards_w_SubjectInfo)

# 'data.frame':	72 obs. of  15 variables:
#   $ Subject.Number       : int  5000 5313 2169 4087 5066 2950 4383 4619 4819 3711 ...
# $ Barcode              : chr  "B007375A" "B010169A" "B000014A" "B010842A" ...
# $ Cohort               : chr  "Cohort 11" "Cohort 13" "Dep Cohort 1" "Dep Cohort 6" ...
# $ Diagnosis            : chr  "Control" "Schiz" "Control" "BP" ...
# $ Age                  : int  77 42 18 51 47 45 40 34 79 59 ...
# $ Gender               : chr  "M" "M" "M" "F" ...
# $ pH                   : num  6.79 6.77 6.97 6.81 6.6 7.05 6.77 6.7 6.75 6.55 ...
# $ AFS                  : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final          : num  16.5 23.7 22 16.8 21 20 26.3 17.2 21 27 ...
# $ Slab.Format          : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number          : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group      : int  1 1 1 1 1 1 2 2 2 2 ...
# $ ID                   : int  1 2 3 4 5 6 7 8 9 10 ...
# $ SampleNumber         : logi  NA NA NA NA NA NA ...
# $ QC_AmplificationIssue: chr  "N" "N" "N" "N" ...


write.csv(Decoder_DA5HTCards_w_SubjectInfo, "Decoder_DA5HTCards_w_SubjectInfo.csv")


Concatenated_DA5HT_wSubjectInfo<-join(Concatenated_DA5HT, Decoder_DA5HTCards_w_SubjectInfo, by="ID", type="left")
str(Concatenated_DA5HT_wSubjectInfo)

# 'data.frame':	7680 obs. of  38 variables:
#   $ Well                       : chr  "A1" "A2" "A3" "A4" ...
# $ ID                         : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Target.Name                : chr  "ADRB1-Hs02330048_s1" "ADRB2-Hs00240532_s1" "COMT-Hs00241349_m1" "DBH-Hs01089840_m1" ...
# $ Amp.Score                  : num  1.25 1.38 1.09 1.31 1.42 ...
# $ Amp.Status                 : chr  "Amp" "Amp" "Amp" "Amp" ...
# $ Task                       : chr  "UNKNOWN" "UNKNOWN" "UNKNOWN" "UNKNOWN" ...
# $ Cq                         : chr  "24.245" "25.087" "24.579" "27.823" ...
# $ Cq.Mean                    : chr  "24.245" "25.087" "24.579" "27.823" ...
# $ Cq.Standard.Deviation      : chr  "-" "-" "-" "-" ...
# $ Quantity                   : chr  "-" "-" "-" "-" ...
# $ Quantity.Mean              : chr  "-" "-" "-" "-" ...
# $ Quantity.Standard.Deviation: chr  "-" "-" "-" "-" ...
# $ Intercept                  : chr  "-" "-" "-" "-" ...
# $ R.Squared                  : chr  "-" "-" "-" "-" ...
# $ Slope                      : chr  "-" "-" "-" "-" ...
# $ Efficiency                 : chr  "-" "-" "-" "-" ...
# $ Auto.Threshold             : chr  "true" "true" "true" "true" ...
# $ Threshold                  : num  0.163 0.239 0.048 0.167 0.24 0.202 0.227 0.194 0.184 0.228 ...
# $ Auto.Baseline              : chr  "true" "true" "true" "true" ...
# $ Baseline.Start             : num  3 3 3 3 3 3 3 3 3 3 ...
# $ Baseline.End               : num  19 20 21 23 27 20 23 27 27 21 ...
# $ Omit                       : chr  "false" "false" "false" "false" ...
# $ X                          : logi  NA NA NA NA NA NA ...
# $ Card                       : chr  "Card 39.eds" "Card 39.eds" "Card 39.eds" "Card 39.eds" ...
# $ Subject.Number             : int  5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 ...
# $ Barcode                    : chr  "B007375A" "B007375A" "B007375A" "B007375A" ...
# $ Cohort                     : chr  "Cohort 11" "Cohort 11" "Cohort 11" "Cohort 11" ...
# $ Diagnosis                  : chr  "Control" "Control" "Control" "Control" ...
# $ Age                        : int  77 77 77 77 77 77 77 77 77 77 ...
# $ Gender                     : chr  "M" "M" "M" "M" ...
# $ pH                         : num  6.79 6.79 6.79 6.79 6.79 6.79 6.79 6.79 6.79 6.79 ...
# $ AFS                        : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final                : num  16.5 16.5 16.5 16.5 16.5 16.5 16.5 16.5 16.5 16.5 ...
# $ Slab.Format                : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number                : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group            : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SampleNumber               : logi  NA NA NA NA NA NA ...
# $ QC_AmplificationIssue      : chr  "N" "N" "N" "N" ...


########################

#Identifying "Housekeeping Genes":

#We need to split the target name column (e.g., "ABAT-Hs00609436_m1") so that the official Gene Symbol is extracted


Temp_delim<- data.frame(do.call('rbind', strsplit(as.character(Concatenated_GabaGlu_wSubjectInfo$Target.Name),'-',fixed=TRUE)))
str(Temp_delim)
# 'data.frame':	13440 obs. of  2 variables:
#   $ X1: Factor w/ 96 levels "18S","ABAT","ACTB",..: 2 4 5 6 7 93 73 10 11 12 ...
# $ X2: Factor w/ 96 levels "Hs00161045_m1",..: 50 62 43 11 45 39 1 87 85 72 ...

Concatenated_GabaGlu_wSubjectInfo$GeneSymbol<-as.character(Temp_delim$X1)
str(Concatenated_GabaGlu_wSubjectInfo)

rm(Temp_delim)

Temp_delim<- data.frame(do.call('rbind', strsplit(as.character(Concatenated_DA5HT_wSubjectInfo$Target.Name),'-',fixed=TRUE)))
str(Temp_delim)

Concatenated_DA5HT_wSubjectInfo$GeneSymbol<-as.character(Temp_delim$X1)
str(Concatenated_GabaGlu_wSubjectInfo)

#We need to use the housekeeping genes list to identify which of the official gene symbols are "housekeeping genes"

Concatenated_GabaGlu_wSubjectInfo$GeneSymbol%in%HousekeepingGenes[,1]

length(Concatenated_GabaGlu_wSubjectInfo$GeneSymbol%in%HousekeepingGenes[,1])
#[1] 13824

sum(Concatenated_GabaGlu_wSubjectInfo$GeneSymbol%in%HousekeepingGenes[,1])
#[1] 1728

Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$GeneSymbol%in%HousekeepingGenes[,1],]
#Excellent - looks correct.

Concatenated_GabaGlu_wSubjectInfo$Housekeeping<-Concatenated_GabaGlu_wSubjectInfo$GeneSymbol%in%HousekeepingGenes[,1]
str(Concatenated_GabaGlu_wSubjectInfo)

Concatenated_GabaGlu_wSubjectInfo$Cq[Concatenated_GabaGlu_wSubjectInfo$Cq=="Undetermined"]<-NA
Concatenated_GabaGlu_wSubjectInfo$Cq.Mean[Concatenated_GabaGlu_wSubjectInfo$Cq.Mean=="Undetermined"]<-NA

Concatenated_GabaGlu_wSubjectInfo$Cq<-as.numeric(Concatenated_GabaGlu_wSubjectInfo$Cq)
Concatenated_GabaGlu_wSubjectInfo$Cq.Mean<-as.numeric(Concatenated_GabaGlu_wSubjectInfo$Cq.Mean)

str(Concatenated_GabaGlu_wSubjectInfo)
write.csv(Concatenated_GabaGlu_wSubjectInfo, "Concatenated_GabaGlu_wSubjectInfo.csv")

Concatenated_DA5HT_wSubjectInfo$Housekeeping<-Concatenated_DA5HT_wSubjectInfo$GeneSymbol%in%HousekeepingGenes[,1]

Concatenated_DA5HT_wSubjectInfo$Cq[Concatenated_DA5HT_wSubjectInfo$Cq=="Undetermined"]<-NA
Concatenated_DA5HT_wSubjectInfo$Cq.Mean[Concatenated_DA5HT_wSubjectInfo$Cq.Mean=="Undetermined"]<-NA

Concatenated_DA5HT_wSubjectInfo$Cq<-as.numeric(Concatenated_DA5HT_wSubjectInfo$Cq)
Concatenated_DA5HT_wSubjectInfo$Cq.Mean<-as.numeric(Concatenated_DA5HT_wSubjectInfo$Cq.Mean)


write.csv(Concatenated_DA5HT_wSubjectInfo, "Concatenated_DA5HT_wSubjectInfo.csv")


#################################################


#QC:

#1. How much do the replicate samples correlate?

#There are 96 measurements per sample on each card.  All cards should follow the same order of targets.
#Sanity Check:
temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="1",]
cbind(temp$GeneSymbol[c(1:96)],temp$GeneSymbol[c(97:192)])
#Yep.
plot(temp$Cq[c(1:96)], temp$Cq[c(97:192)])
#Ooh - nice correlation.
#Although it has one super high leverage point-18S has very low Cq in both samples (highly expressed) 

#Percentage difference between the cards for each target:
(abs(temp$Cq[c(1:96)]-temp$Cq[c(97:192)])/apply(cbind(temp$Cq[c(1:96)],temp$Cq[c(97:192)]), 1, mean, na.rm=T))*100

#Let's loop it:

length(names(table(Concatenated_GabaGlu_wSubjectInfo$ID)))
#[1] 72


ReplicateCor_ForEachSubjectID_GabaGlu<-matrix(0, 72, 1)
row.names(ReplicateCor_ForEachSubjectID_GabaGlu)<-names(table(Concatenated_GabaGlu_wSubjectInfo$ID))

Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates<-rep(0, length(Concatenated_GabaGlu_wSubjectInfo[,1]))

for (i in c(1:72)){
  
  print(names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i])
  
  temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID==names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i],] 

  ReplicateCor_ForEachSubjectID_GabaGlu[i]<-cor(x=temp$Cq[c(1:96)], y=temp$Cq[c(97:192)], use="pairwise.complete.obs")
  
  PercentageDiff_BetweenReplicates<-(abs(temp$Cq[c(1:96)]-temp$Cq[c(97:192)])/apply(cbind(temp$Cq[c(1:96)],temp$Cq[c(97:192)]), 1, mean, na.rm=T))*100
  print(length(PercentageDiff_BetweenReplicates))
  print(length(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_GabaGlu_wSubjectInfo$ID==names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i]] ))
  Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_GabaGlu_wSubjectInfo$ID==names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i]] <-rep(PercentageDiff_BetweenReplicates, 2)
  
  rm(temp)
  rm(PercentageDiff_BetweenReplicates)
  
}

#Original notes - this problem is now fixed:
#Ah ha - I found a problem. "Card 22.eds" isn't in my PCR data folder. Interesting. That means two things:
##1) My decoder is wrong regarding which subject IDs are found on which cards (a minor point, since I didn't end up using it, just curious).
##2) I'm missing a file.  I double-checked the server and it isn't there either - so not a download issue on my end. E-mailed Adriana.
##3) She re-uploaded the the missing file, and I re-ran the analysis on 08-25-2020 


#27 just says NA for the correlation, but has the correct number of genes.  
i<-19
print(names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i])
[1] "27"
temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID==names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i],]
View(temp)
#Ah - apparently the second run of the data (rows temp$Cq[c(97:192)]) includes exactly one measurement (18S - the most highly expressed gene) and everything else is NA. Ouch. 
#Which is weird, because this wasn't one of the cards identified by Adriana as having a QC issue.
Decoder_GabaGluCards 
#Subject.Number ID      Folder SampleNumber QC_AmplificationIssue
#52            4301 27 Card 14.eds            2                     Y

#Ah - actually it was identified as having a QC issue originally but that replicate-specific info was lost during the joining process (because we joined by subject ID)


#Let's snoop at it:
print(names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i])

write.csv(ReplicateCor_ForEachSubjectID_GabaGlu, "ReplicateCor_ForEachSubjectID_GabaGlu.csv")

#...and 38 says cor=1... which is impossible unless the data are identical. So there's a bug here.
#It's not a problem with sample id formatting - there are 72 IDs and I don't see doubles for samples 27, 41-44 with slightly different formatted IDs.
Decoder_GabaGluCards
#Subject.Number ID      Folder SampleNumber QC_AmplificationIssue
#69            4359 38 Card 17.eds            3                     Y
#This is also one of the samples that had a QC issue.


#In the meantime I'm going to snoop at the data that I do have:

boxplot(ReplicateCor_ForEachSubjectID_GabaGlu[,1])
#Pretty tight. Whiskers & IQ all fall above 0.98 (maybe more like 0.985?). There are 10 "outlier" datapoints outside of that range. Let's see who they are:
ReplicateCor_ForEachSubjectID_GabaGlu[ReplicateCor_ForEachSubjectID_GabaGlu[,1]<0.984,]
# 1        13        15        16        25      <NA>        37        49        61        69         7 
# 0.8651958 0.9583229 0.9795163 0.9715362 0.9652474        NA 0.9737637 0.9398720 0.9773817 0.9835464 0.9802135 

#Why is there a name NA?
rownames(ReplicateCor_ForEachSubjectID_GabaGlu)#not NA here
#I think it must be #27 though, since that is the only sample with no duplicate.

#Other notes:
#Sample ID #1 has a particularly low replicate sample correlation (0.86)
#49 and 13 also have notably low replicate sample correlations.

#It would be helpful to interpret this data in the context of the correlations across subjects (vs. within subjects). Our data frame is currently not well set up for that 

#If I could rearrange it, I could reuse a lot of my code from microarray/RNA-Seq analyses.

dim(Concatenated_GabaGlu_wSubjectInfo)
#[1] 13824    41

13824/96
#[1] 144
144/2
#[1] 72

Concatenated_GabaGlu_wSubjectInfo$ID[c(1:96)]
Concatenated_GabaGlu_wSubjectInfo$ID[c(97:(96+96))]

table(Concatenated_GabaGlu_wSubjectInfo$ID[seq(1, 13824, 96)])
#looks good: two samples per subject

Concatenated_GabaGlu_wSubjectInfo$ID[seq(1, 13824, 96)]


GabaGlu_Cq_AllSubjects<-matrix(0,96,144)
colnames(GabaGlu_Cq_AllSubjects)<-Concatenated_GabaGlu_wSubjectInfo$ID[seq(1, 13824, 96)]
colnames(GabaGlu_Cq_AllSubjects)
#If we use the subject IDs, the columnames sometimes repeat - I thought R would be unhappy about it, but I don't see a .2 added to the names. Interesting.

FirstGeneForSample_iterative<-seq(1, 13824, 96)
LastGeneForSample_iterative<-seq(96, 13824, 96)

i<-1
#Sanity check:
Concatenated_GabaGlu_wSubjectInfo$ID[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]
row.names(GabaGlu_Cq_AllSubjects)<-Concatenated_GabaGlu_wSubjectInfo$GeneSymbol[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]


for(i in c(1:144)){
  GabaGlu_Cq_AllSubjects[,i]<-Concatenated_GabaGlu_wSubjectInfo$Cq[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]
  
}

head(GabaGlu_Cq_AllSubjects)
#yesss....
#colnames 27(v2) and 38 (v2) are mostly NA data.
colnames(GabaGlu_Cq_AllSubjects)[23]
colnames(GabaGlu_Cq_AllSubjects)[50]
#Note: this will need to be updated once we sort out the issue with the missing file***********

SubjectInfo_OrderedForGabaGluCqMatrix<-Concatenated_GabaGlu_wSubjectInfo[FirstGeneForSample_iterative, c(2, 24:38)]
head(SubjectInfo_OrderedForGabaGluCqMatrix)

GabaGlu_Cq_AllSubjects_QCed<-GabaGlu_Cq_AllSubjects[,-c(23,50)]
SubjectInfo_OrderedForGabaGluCqMatrix_QCed<-SubjectInfo_OrderedForGabaGluCqMatrix[-c(23,50), ]
#Sanity check:
colnames(GabaGlu_Cq_AllSubjects_QCed)
SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID
#Yep, bad samples 27(v2) and 38 (v2) now gone & ID orders match in the two dataframes.


GabaGlu_SubjectSubjectCorMatrix<-cor(GabaGlu_Cq_AllSubjects_QCed, use="pairwise.complete.obs")

pdf("GabaGlu_SubjectSubjectCorMatrix.pdf", height=14, width=14)
heatmap(GabaGlu_SubjectSubjectCorMatrix)
dev.off()
#Interesting - some of the replicates most closely resemble each other, but certainly not all (maybe 20%?). There are correlation blocks (batch effects??). So the data is definitely noisy.

#To do:
#Try spearman
#Make an example plot
#Run PCA
#Examine influence of other variables (e.g., chip)
#Try other dataset


#Example plots:
plot(GabaGlu_Cq_AllSubjects_QCed[,1]~GabaGlu_Cq_AllSubjects_QCed[,2])
plot(GabaGlu_Cq_AllSubjects_QCed[,4]~GabaGlu_Cq_AllSubjects_QCed[,5])
#There is definitely one super high leverage point with Cq less than 10 (18S - housekeeping gene with super high expression) - I might try redoing the correlation analysis with that point removed and/or z-score the data first so that the relative values across subjects are emphasized instead of the the higher leverage points.

GabaGlu_Cq_AllSubjects_Scaled_QCed<-t(scale(t(GabaGlu_Cq_AllSubjects_QCed), center=T, scale=T))
head(GabaGlu_Cq_AllSubjects_Scaled_QCed)
#Sanity check:
apply(GabaGlu_Cq_AllSubjects_Scaled_QCed, 1, mean)
apply(GabaGlu_Cq_AllSubjects_Scaled_QCed, 1, sd)
#Looks good.

GabaGlu_SubjectSubjectCorMatrix_Scaled_QCed<-cor(GabaGlu_Cq_AllSubjects_Scaled_QCed, use="pairwise.complete.obs")
GabaGlu_SubjectSubjectCorMatrix_Scaled_QCed

pdf("GabaGlu_SubjectSubjectCorMatrix_Scaled_QCed.pdf", height=14, width=14)
heatmap(GabaGlu_SubjectSubjectCorMatrix_Scaled_QCed)
dev.off()
#The correlation blocks are even more obvious now- let's see what variables they relate to

# GabaGlu_SubjectSubjectCorMatrix_Spearman<-cor(GabaGlu_Cq_AllSubjects_QCed, use="pairwise.complete.obs", method="spearman")

# pdf("GabaGlu_SubjectSubjectCorMatrix_Spearman.pdf", height=14, width=14)
# heatmap(GabaGlu_SubjectSubjectCorMatrix_Spearman)
# dev.off()
# #The replicates are less obvious using Spearman


#pca_GabaGlu<-prcomp(na.omit(t(GabaGlu_Cq_AllSubjects_Scaled_QCed)))
#I originally used this code because I was getting errors due to NA values and discovered that 100 of my subjects were omitted - ouch!  Let's try a different version.

sum(is.na(GabaGlu_Cq_AllSubjects_Scaled_QCed))
#[1] 133
length(GabaGlu_Cq_AllSubjects_Scaled_QCed)
#[1] 13632
#Percent NA:
(sum(is.na(GabaGlu_Cq_AllSubjects_Scaled_QCed))/length(GabaGlu_Cq_AllSubjects_Scaled_QCed))*100
#[1] 0.9756455
#Less than 1% NA measurements.

#Since this is scale data, the mean value is 0 for all genes. Let's replace NAs with that:
GabaGlu_Cq_AllSubjects_Scaled_QCed[is.na(GabaGlu_Cq_AllSubjects_Scaled_QCed)]<-0

pdf("Boxplot_CqZscore_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for all genes", xlab="Sample ID", ylab="Cq Z-score")
dev.off()
#Huge differences in relative Cq across samples. Seems to partially correlate with card. 
#So we definitely need to do something to normalize for overall signal levels on the cards (which is why PCR normally normalizes by housekeeping genes)
#... but that won't help the lack of correlation between replicate samples, since it only corrects for overall Cq differences across samples.

#For comparison with later -DeltaCq plot:
pdf("Boxplot_CqZscore_perSample_AllGenes_ylim10.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for all genes", xlab="Sample ID", ylab="Cq Z-score", ylim=c(-10,10))
dev.off()

#Let's look at the "housekeeping genes" specifically:

cbind(row.names(GabaGlu_Cq_AllSubjects_Scaled_QCed), Concatenated_GabaGlu_wSubjectInfo$Housekeeping[c(1:96)])
#Looks like the housekeeping genes are in well 11 (18S) and the last 10 wells for each sample (86-96)

pdf("Boxplot_CqZscore_perSample_HousekeepingGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_Scaled_QCed)[c(11, 86:96),], cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for all genes", xlab="Sample ID", ylab="Cq Z-score")
dev.off()
#Yes, the housekeeping genes are tightly variable across samples. They definitely seem to vary with card, but also vary by sample in a manner that does not necessarily correlate across replicates.

GabaGlu_TrimmedMeanHousekeepingZScore<-apply(GabaGlu_Cq_AllSubjects_Scaled_QCed[c(11,86:96),], 2, function(y) mean(y, trim=0.2))
#Side note: This is just me screwing around.
hist(GabaGlu_TrimmedMeanHousekeepingZScore)


#Trying PCA again:
pca_GabaGlu<-prcomp(t(GabaGlu_Cq_AllSubjects_Scaled_QCed))
tmp<-pca_GabaGlu$x[,1:10]
dim(tmp)
#[1] 142  10
rownames(tmp)<-colnames(GabaGlu_Cq_AllSubjects_Scaled_QCed)
write.csv(tmp, "PCA_GabaGlu.csv")

tmp<-pca_GabaGlu$rotation[,1:10]
write.csv(tmp, "pca_GabaGlu_Eigenvectors.csv")

png("PCA_ScreePlot_GabaGlu.png")
plot(summary(pca_GabaGlu)$importance[2,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 is a huge amount of the variance (>60%)

png("PCA_ScreePlot2_GabaGlu.png")
plot(summary(pca_GabaGlu)$importance[3,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()


png("PC1vsPC2_GabaGlu_byDiagnosis.png")
plot(pca_GabaGlu$x[,1]~pca_GabaGlu$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2", col=as.factor(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis))
dev.off()
#Not diagnosis.
#Also not driven by a single major outlier.

png("PC1_GabaGlu_byCohort.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Cohort, xlab="Cohort", ylab="PC1")
dev.off()
#Varies a little bit by cohort, but not strikingly so.

pdf("PC1_GabaGlu_byCard.pdf", width=10, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card, xlab="Cohort", ylab="PC1", cex.lab=0.3, las=3)
dev.off()
#Yeah, looks like there might be a fair amount of variation based on card, but it isn't everything.
summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card))
# Multiple R-squared:  0.3626,	Adjusted R-squared:  0.1521 
# F-statistic: 1.723 on 35 and 106 DF,  p-value: 0.01813
#Cards that are particularly odd:
#10, 18, 20, 22, 34, 36, 6 (especially 18 &  6)
#Since the same samples were in adjacent card numbers (e.g., 1, 2 - although note that these are not adjacent in the plot), it is easy to visually see that these card differences aren't due to differences in sample make-up - it is actually the cards themselves. 
#Since there are only 4 samples per card, we probably can't explicitly correct for this (maybe with a multilevel model?)
#It should be partially helped by averaging the replicates if we go that route.

#Out of curiousity, let's see if correcting for it improves clustering of replicates:

pdf("PC1_GabaGlu_byID.pdf", width=10, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID, ylim=c(-15, 45))
dev.off()

str(summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card)))

summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card))$residuals

pdf("PC1_GabaGlu_byID_CardEffectRemoved.pdf", width=10, height=5)
boxplot(summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card))$residuals~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID, ylim=c(-15, 45))
dev.off()
#Yeah, that definitely tightens things up.
#Next decision - just average by card or explicitly include it in a multilevel model?

pdf("PC1_GabaGlu_byTrimmedMeanHousekeepingZScore.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~GabaGlu_TrimmedMeanHousekeepingZScore)
dev.off()
#Yep, very interesting - the main source of variation in the data is just some overall effect on Cq, which is easily controlled using some sort of mean Cq for the Housekeeping Genes. That's strange though, since the PC scores are based on the correlation across samples (which shouldn't be affected by overall Cq, unless overall Cq is having some sort of non-linear effect on the signal for some of the genes)


png("PC1_GabaGlu_byDissection.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Dissecton.Group, xlab="Dissecton.Group", ylab="PC1")
dev.off()
#Some variation related to dissection, but not strikingly so.

png("PC1_GabaGlu_byGender.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Gender, xlab="Gender", ylab="PC1")
dev.off()
#Not Gender

png("PC1_GabaGlu_byAge.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Age, xlab="Age", ylab="PC1")
dev.off()
#PC1 is not age


#Some of the pH's are recorded as 0 - that's problematic.
#cbind(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Subject.Number, SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH)
#Subject #4463 is the one missing pH info.
#I double checked with the 2019 version of the Pritzker database that I have, and it is missing there too. I'll need to double check on the server, and then probably contact UCIrvine.
#For now, let's replace with NA
#Since we are going to need to re-run this analysis anyway, I went back to the original .csv document and replaced with NA, but for now I'll just fix this in the current Subject Info matrix to peek at our data.
#SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH[SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH==0]<-NA

png("PC1_GabaGlu_byPH.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH, xlab="pH", ylab="PC1")
dev.off()
#PC1 is not pH, although it may slightly correlate

colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed)

png("PC1_GabaGlu_byHoursFinal.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Hours.Final, xlab="Hours.Final", ylab="PC1")
dev.off()
#PC1 is not PMI, maybe a slight correlation

png("PC1_GabaGlu_byID.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID, xlab="ID", ylab="PC1")
dev.off()
#PC1 is not ID, maybe a slight correlation

png("PC1_GabaGlu_bySubjectNumber.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Subject.Number, xlab="Subject.Number", ylab="PC1")
dev.off()
#No linear relationship with subject number


#Checking to see if the housekeeping genes vary with diagnosis (note: this hasn't been averaged by chip yet)

#Making sure Diagnosis has the correct reference group:
levels(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
#SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis<-as.factor(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
#SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis<-relevel(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis, ref="Control")
#levels(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
# 
# for(i in c(86:96)){
#   print(row.names(GabaGlu_Cq_AllSubjects_QCed)[i])
# print(summary.lm(lm(GabaGlu_Cq_AllSubjects_QCed[i,]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)))
# }
# #Varies with Schiz: GUSB (p=0.0256),  IP08 (p=0.0207), TBP (p=0.0417), TFRC (p=0.0328)
# #Note: Card matters a lot.
# 
# #Effect of diagnosis on housekeeping genes, after controlling for card:
# for(i in c(86:96)){
#   print(row.names(GabaGlu_Cq_AllSubjects_QCed)[i])
#   print(summary.lm(lm(GabaGlu_Cq_AllSubjects_QCed[i,]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis+SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card)))
# }
# #Varies with Schiz: GUSB (p=0.0202),  IP08 (p=0.0354), TBP (p=0.0647), TFRC (p=0.044)
# #Varies with BP: HMBS(p=0.04324)
# 
# #checking on the usual suspects:
# for(i in c(86:96)){
#   print(row.names(GabaGlu_Cq_AllSubjects_QCed)[i])
#   print(summary.lm(lm(GabaGlu_Cq_AllSubjects_QCed[i,]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis+SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Age+SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Gender+SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH+SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Hours.Final+SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card)))
# }
# #Age seems to be particularly important as a co-variate. Even in this large of a model, Age is significantly related to housekeeping gene expression for 7/10 housekeeping genes
# 
# summary.lm(lm(SubjectInfo$Age~SubjectInfo$Diagnosis))
# #Call:
# # lm(formula = SubjectInfo$Age ~ SubjectInfo$Diagnosis)
# # 
# # Residuals:
# #   Min      1Q  Median      3Q     Max 
# # -33.370  -8.010  -0.519   6.806  27.630 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)                  51.370      2.567  20.008   <2e-16 ***
# #   SubjectInfo$DiagnosisBP      -5.704      3.882  -1.469    0.146    
# # SubjectInfo$DiagnosisSchiz   -4.912      3.743  -1.312    0.194    
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Residual standard error: 13.34 on 69 degrees of freedom
# # Multiple R-squared:  0.03744,	Adjusted R-squared:  0.009544 
# # F-statistic: 1.342 on 2 and 69 DF,  p-value: 0.268
# 
# #Age tends to be younger in the diagnosis subjects by about 5 years, although not significantly so.
# 
# summary.lm(lm(SubjectInfo$pH~SubjectInfo$Diagnosis))
# # Call:
# #   lm(formula = SubjectInfo$pH ~ SubjectInfo$Diagnosis)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.24381 -0.12313 -0.01381  0.11735  0.35619 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)                 6.831481   0.029050 235.163   <2e-16 ***
# #   SubjectInfo$DiagnosisBP    -0.067672   0.043920  -1.541    0.128    
# # SubjectInfo$DiagnosisSchiz -0.006699   0.042832  -0.156    0.876    
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Residual standard error: 0.1509 on 68 degrees of freedom
# # (1 observation deleted due to missingness)
# # Multiple R-squared:  0.03863,	Adjusted R-squared:  0.01035 
# # F-statistic: 1.366 on 2 and 68 DF,  p-value: 0.262
# 
# #pH may be a little lower in the BP subjects, but not significantly so.
# 
# 
# summary.lm(lm(SubjectInfo$Hours.Final~SubjectInfo$Diagnosis))
# summary.lm(lm(SubjectInfo$Hours.Final~SubjectInfo$Diagnosis))
# 
# # Call:
# #   lm(formula = SubjectInfo$Hours.Final ~ SubjectInfo$Diagnosis)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -13.3143  -5.5717  -0.2229   4.5978  23.8062 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)                 22.7315     1.3698  16.595   <2e-16 ***
# #   SubjectInfo$DiagnosisBP     -0.4172     2.0709  -0.201    0.841    
# # SubjectInfo$DiagnosisSchiz  -2.2877     1.9968  -1.146    0.256    
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Residual standard error: 7.117 on 69 degrees of freedom
# # Multiple R-squared:  0.02045,	Adjusted R-squared:  -0.007945 
# # F-statistic: 0.7202 on 2 and 69 DF,  p-value: 0.4903
# 
# 
# table(SubjectInfo$Gender, SubjectInfo$Diagnosis)
# # Control BP Schiz
# # F       1  6     1
# # M      26 15    23
# fisher.test(SubjectInfo$Gender, SubjectInfo$Diagnosis)
# # Fisher's Exact Test for Count Data
# # 
# # data:  SubjectInfo$Gender and SubjectInfo$Diagnosis
# # p-value = 0.01363
# # alternative hypothesis: two.sided

#Alright - we'll definitely have to control for that one. 

library(nlme)
i<-1
Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
str(Temp)
summary(lme(y~Diagnosis+Age+Gender+pH+Hours.Final, random=~Card|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim')))
#This code crashes R. I think I must be including Card incorrectly in the model.

summary(lme(y~Diagnosis+Age+Gender+pH+Hours.Final, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim')))

summary(lme(y~Diagnosis+Age+Gender+pH+Hours.Final+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim')))
#That works, but I would really prefer to have Card as a random effect (since Card is just a batch effect)

#Ah - here we go - I did have the nested terms specified wrong: https://biostatmatt.com/archives/2718
#With proper specification:
summary(lme(y~Diagnosis+Age+Gender+pH+Hours.Final, random=~1|Card/ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim')))
#Or maybe that isn't quite right?  There are two samples for each ID found on two different cards. So technically they are independent (crossed) random variables, sort-of...
#https://errickson.net/stats-notes/vizrandomeffects.html
#This seems to suggest that it is difficult to specify crossed random effects in nlme and that I should switch to lme4 instead. Sigh:
#https://biostatmatt.com/archives/2718


library(lme4)
lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | Card) + (1 | ID), data = Temp)
#Very similar to what is given above using nlme... 
summary(lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | Card) + (1 | ID), data = Temp))
#except that there are no p-values. I remember this being a problem last time I used lme4 too, and that the work around for it was weird. Sigh.
#https://rdrr.io/cran/lme4/man/pvalues.html

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | Card) + (1 | ID), data = Temp)
summary(Model)

# Linear mixed model fit by REML ['lmerMod']
# Formula: y ~ Diagnosis + Age + Gender + pH + Hours.Final + (1 | Card) +      (1 | ID)
# Data: Temp
# 
# REML criterion at convergence: 196.8
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.8282 -0.5157 -0.1561  0.4931  2.5193 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# ID       (Intercept) 0.06999  0.2646  
# Card     (Intercept) 0.02307  0.1519  
# Residual             0.12508  0.3537  
# Number of obs: 140, groups:  ID, 71; Card, 36
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)    23.767696   2.129299  11.162
# DiagnosisBP    -0.216930   0.118327  -1.833
# DiagnosisSchiz -0.025565   0.106535  -0.240
# Age            -0.006668   0.003493  -1.909
# GenderM         0.035745   0.153541   0.233
# pH             -0.357061   0.307838  -1.160
# Hours.Final     0.003967   0.006544   0.606
# 
# Correlation of Fixed Effects:
#             (Intr) DgnsBP DgnssS Age    GendrM pH    
# DiagnosisBP -0.235                                   
# DiagnssSchz -0.064  0.422                            
# Age         -0.199  0.181  0.117                     
# GenderM      0.006  0.305 -0.026 -0.012              
# pH          -0.993  0.186  0.026  0.124 -0.062       
# Hours.Final -0.030 -0.070  0.119 -0.098 -0.174 -0.019

car::Anova(Model)

# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# Diagnosis   3.7085  2    0.15657  
# Age         3.6439  1    0.05627 .
# Gender      0.0542  1    0.81591  
# pH          1.3454  1    0.24609  
# Hours.Final 0.3675  1    0.54437  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 124.5950  1    < 2e-16 ***
#   Diagnosis     3.7085  2    0.15657    
# Age           3.6439  1    0.05627 .  
# Gender        0.0542  1    0.81591    
# pH            1.3454  1    0.24609    
# Hours.Final   0.3675  1    0.54437    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting - we're getting exactly the same answer with Type II and Type III. Huh.

#https://www.researchgate.net/post/Problems_performing_ANOVA_and_Mixed_Linear_model_on_R
#https://www.researchgate.net/post/REML_FALSE_versus_REML_TRUE_lme4_package_in_R-any_thoughts
#Apparently REML (which is the standard for evaluating multilevel models) shouldn't be used in combination with anova, but anova() function may automatically correct it (make it ML)
#Let's check:

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | Card) + (1 | ID), data = Temp, REML=F)
summary(Model)

# Linear mixed model fit by maximum likelihood  ['lmerMod']
# Formula: y ~ Diagnosis + Age + Gender + pH + Hours.Final + (1 | Card) +      (1 | ID)
# Data: Temp
# 
# AIC      BIC   logLik deviance df.resid 
# 186.6    216.0    -83.3    166.6      130 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.9471 -0.4970 -0.1629  0.5083  2.6107 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# ID       (Intercept) 0.05615  0.2370  
# Card     (Intercept) 0.02041  0.1429  
# Residual             0.12679  0.3561  
# Number of obs: 140, groups:  ID, 71; Card, 36
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)    23.780109   2.021971  11.761
# DiagnosisBP    -0.217862   0.112429  -1.938
# DiagnosisSchiz -0.025262   0.101188  -0.250
# Age            -0.006689   0.003318  -2.016
# GenderM         0.035288   0.145776   0.242
# pH             -0.358705   0.292297  -1.227
# Hours.Final     0.003982   0.006212   0.641
# 
# Correlation of Fixed Effects:
#   (Intr) DgnsBP DgnssS Age    GendrM pH    
# DiagnosisBP -0.236                                   
# DiagnssSchz -0.064  0.422                            
# Age         -0.199  0.182  0.117                     
# GenderM      0.005  0.305 -0.025 -0.011              
# pH          -0.993  0.186  0.025  0.124 -0.062       
# Hours.Final -0.030 -0.069  0.120 -0.098 -0.174 -0.019

car::Anova(Model, type="III")

# Linear mixed model fit by maximum likelihood  ['lmerMod']
# Formula: y ~ Diagnosis + Age + Gender + pH + Hours.Final + (1 | Card) +      (1 | ID)
# Data: Temp
# 
# AIC      BIC   logLik deviance df.resid 
# 186.6    216.0    -83.3    166.6      130 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.9471 -0.4970 -0.1629  0.5083  2.6107 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# ID       (Intercept) 0.05615  0.2370  
# Card     (Intercept) 0.02041  0.1429  
# Residual             0.12679  0.3561  
# Number of obs: 140, groups:  ID, 71; Card, 36
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)    23.780109   2.021971  11.761
# DiagnosisBP    -0.217862   0.112429  -1.938
# DiagnosisSchiz -0.025262   0.101188  -0.250
# Age            -0.006689   0.003318  -2.016
# GenderM         0.035288   0.145776   0.242
# pH             -0.358705   0.292297  -1.227
# Hours.Final     0.003982   0.006212   0.641
# 
# Correlation of Fixed Effects:
#   (Intr) DgnsBP DgnssS Age    GendrM pH    
# DiagnosisBP -0.236                                   
# DiagnssSchz -0.064  0.422                            
# Age         -0.199  0.182  0.117                     
# GenderM      0.005  0.305 -0.025 -0.011              
# pH          -0.993  0.186  0.025  0.124 -0.062       
# Hours.Final -0.030 -0.069  0.120 -0.098 -0.174 -0.019
# > car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 138.3177  1    < 2e-16 ***
#   Diagnosis     4.1471  2    0.12574    
# Age           4.0648  1    0.04379 *  
#   Gender        0.0586  1    0.80872    
# pH            1.5060  1    0.21975    
# Hours.Final   0.4111  1    0.52143    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#############

#Side note: NA omit means that the subject with missing pH may not be included in the analysis - and in general, the results may change once they are included. So this analysis is preliminary.
length(Temp[,1])
#[1] 142
#whereas above the Number of obs: 140 
#Maybe we should replace with average pH in the meantime?

SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH[is.na(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH)]<-mean(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH, na.rm=T)
str(SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Dissecton.Group<-as.factor(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Dissecton.Group)
str(SubjectInfo_OrderedForGabaGluCqMatrix_QCed)

Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
str(Temp)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | Card) + (1 | ID), data = Temp, REML=F)
summary(Model)
#Now we have all 142 observations. :)

#Alright - let's snoop at those housekeeping genes again:
#The usual suspects:


for(i in c(11,86:96)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspects.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
    cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspects.txt", sep="\n", append=TRUE)
    close(OutputtingStats)
    rm(stats_output)
  rm(Temp)
  rm(Model)
}

#Age sometimes matters, Hours Final sometimes matters, Diagnosis does not significantly matter (trend for GUSB)

#The usual suspects + some batch related variables:



table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Cohort, SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
#                 Control BP Schiz
# Cohort 11           12  2     6
# Cohort 12           18  8     4
# Cohort 13           10  0     4
# Cohort 7             0 11    11
# Cohort 8             0  4     8
# Dep Cohort 1         4  4     0
# Dep Cohort 5         2  8     0
# Dep Cohort 6         0  4     4
# Schiz Cohort 2       8  0    10

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Dissecton.Group, SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Cohort, SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
#       Control BP Schiz
# 1        6  2     4
# 2        4  4     2
# 3        4  4     4
# 4        4  4     4
# 5        4  4     3
# 6        4  2     6
# 7        6  1     4
# 8        4  4     4
# 9        4  4     4
# 10       4  4     4
# 11       4  4     4
# 12       6  4     4

#Dissection group is reasonably well balanced across diagnosis (go me!), cohort is not.

for(i in c(11, 86:96)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Cohort+ (1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndBatchFixed.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndBatchFixed.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Dissection group and cohort matter for many of these genes. After controlling for them, diagnosis matters for one gene (HPRT1)

#Side note: Dissection batch and Cohort could also be modeled as random effects (instead of fixed)
#If so, ID is nested in Cohort and in Dissection Group
table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Cohort, SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Dissecton.Group)
table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Cohort, SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card)
#               1 2 3 4 5 6 7 8 9 10 11 12
# Cohort 11      2 2 2 0 2 2 2 2 2  0  2  2
# Cohort 12      2 2 2 2 4 4 2 4 0  2  4  2
# Cohort 13      2 0 0 2 0 0 2 2 4  0  0  2
# Cohort 7       0 2 2 4 3 2 1 2 0  4  0  2
# Cohort 8       0 2 2 0 0 0 2 0 2  2  0  2
# Dep Cohort 1   2 0 0 0 0 0 0 2 2  0  0  2
# Dep Cohort 5   0 2 2 2 0 0 0 0 0  2  2  0
# Dep Cohort 6   2 0 0 2 0 0 0 0 2  0  0  2
# Schiz Cohort 2 2 0 2 0 2 4 2 0 0  2  4  0

#So how do we set that up?
#https://biostatmatt.com/archives/2718

SubjectInfo_OrderedForGabaGluCqMatrix_QCed$CohortDissection<-with(SubjectInfo_OrderedForGabaGluCqMatrix_QCed, interaction(Cohort,Dissecton.Group))
table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$CohortDissection)
#Most of these combinations only include 2 samples (i.e., one sample ID) - not very useful.


Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | Card) + (1 | Cohort/ID)+(1|Dissecton.Group/ID), data = Temp, REML=F)
summary(Model)
#I don't think that is correct.

#https://errickson.net/stats-notes/vizrandomeffects.html
#We can include nested random effects using the cross effects syntax. 


Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1|Dissecton.Group)+(1|Cohort)+(1 | Card) + (1 | ID), data = Temp, REML=F)
summary(Model)
car::Anova(Model, type="III")
#The Stdev assigned to Cohort is 0 - something seems wrong still. Maybe because Cohort*Dissection Group is almost the same thing as ID? But we're treating the effects as independent, so I'm not sure why that would arise as an issue.

#...but since cohort and diagnosis are correlated, perhaps we should include cohort as a fixed effect?
Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
summary(Model)
car::Anova(Model, type="III")
#hmmm.... but since diagnosis and cohort are not completely independent, do we need to worry about how they are centered?  Normally linear models are better with unbalanced design, right? (i.e., so that the effect of cohort is estimated within control subjects)


#This version with cohort as a fixed effect definitely seems to more closely parallel the results with cohort and dissection as fixed effects
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Cohort+ (1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")


for(i in c(11,86:96)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndBatchRandom.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndBatchRandom.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#HPRT1 is the only housekeeping gene with a significant diagnosis effect. Several genes have significant effects of cohort, age, or PMI

#So let's try controlling for average housekeeping gene Cq to clean up the data:

sum(is.na(GabaGlu_Cq_AllSubjects_QCed[c(11,86:96),]))
#[1] 0
#Nice - I don't have to decide what to do with the NAs...
hist(GabaGlu_Cq_AllSubjects_QCed[c(11,86:96),])
#Range Cq=5-29, somewhat normally distributed except for 18S (which has really low Cq)
row.names(GabaGlu_Cq_AllSubjects_QCed[c(11,86:96),])

#I'm going to try to use the trimmed mean for the reference gene signal so that 18S isn't quite as influential:

GabaGlu_MeanHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed[c(11, 86:96),], 2, function(y) mean(y))

GabaGlu_TrimmedMeanHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed[c(11, 86:96),], 2, function(y) mean(y, trim=0.2))

hist(GabaGlu_TrimmedMeanHousekeeping)

#They're almost identical - so I guess either one would probably be fine.
plot(GabaGlu_TrimmedMeanHousekeeping~GabaGlu_MeanHousekeeping)

GabaGlu_TrimmedMeanNonHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed[-c(11, 86:96),], 2, function(y) mean(y, trim=0.2, na.rm=T))

GabaGlu_MeanNonHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed[-c(11, 86:96),], 2, function(y) mean(y, na.rm=T))

pdf("GabaGlu_MeanNonHousekeepingVsMeanHousekeeping.pdf", height=4, width=4)
plot(GabaGlu_MeanNonHousekeeping~GabaGlu_MeanHousekeeping)
dev.off()
cor(GabaGlu_MeanNonHousekeeping,GabaGlu_MeanHousekeeping)
#[1] 0.8654978

pdf("GabaGlu_TrimmedMeanNonHousekeepingVsTrimmedMeanHousekeeping.pdf", height=4, width=4)
plot(GabaGlu_TrimmedMeanNonHousekeeping~GabaGlu_TrimmedMeanHousekeeping)
dev.off()
cor(GabaGlu_TrimmedMeanNonHousekeeping,GabaGlu_TrimmedMeanHousekeeping)
#[1] 0.870455

#Mean housekeeping and non-housekeeping are really strongly correlated.


#Calculating -Delta Cq (Target-Housekeeping):

boxplot((GabaGlu_Cq_AllSubjects_QCed[1,]*-1)~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
boxplot(((GabaGlu_Cq_AllSubjects_QCed[1,]-GabaGlu_MeanHousekeeping)*-1)~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
#There are some extreme outliers there.
boxplot(((GabaGlu_Cq_AllSubjects_QCed[1,]-GabaGlu_TrimmedMeanHousekeeping)*-1)~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)

boxplot((GabaGlu_Cq_AllSubjects_QCed[2,]*-1)~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
boxplot(((GabaGlu_Cq_AllSubjects_QCed[2,]-GabaGlu_MeanHousekeeping)*-1)~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
boxplot(((GabaGlu_Cq_AllSubjects_QCed[2,]-GabaGlu_TrimmedMeanHousekeeping)*-1)~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)

#Visually Neither mean nor trimmed mean seems to perform better.


head(GabaGlu_Cq_AllSubjects_QCed)

GabaGlu_NegDeltaCq_AllSubjects_QCed<-apply(GabaGlu_Cq_AllSubjects_QCed, 1, function(y) (y-GabaGlu_MeanHousekeeping)*-1)
str(GabaGlu_NegDeltaCq_AllSubjects_QCed)
# num [1:142, 1:96] 0.374 0.152 0.276 0.386 0.404 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:142] "1" "2" "3" "4" ...
# ..$ : chr [1:96] "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
head(GabaGlu_NegDeltaCq_AllSubjects_QCed)
#The apply function transposed things. Huh. I wonder why.

GabaGlu_TrimmedMeanHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed[c(11, 86:96),], 2, function(y) mean(y, trim=0.2))

GabaGlu_NegDeltaCqTrimmed_AllSubjects_QCed<-apply(GabaGlu_Cq_AllSubjects_QCed, 1, function(y) (y-GabaGlu_TrimmedMeanHousekeeping)*-1)
str(GabaGlu_NegDeltaCqTrimmed_AllSubjects_QCed)


sum(is.na(GabaGlu_Cq_AllSubjects_QCed))
#[1] 133

pdf("Boxplot_Cq_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq per sample for all genes", xlab="Sample ID", ylab="Cq")
dev.off()

sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_QCed))

pdf("Boxplot_NegDeltaCq_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(t(GabaGlu_NegDeltaCq_AllSubjects_QCed)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq per sample for all genes", xlab="Sample ID", ylab="-DeltaCq")
dev.off()

pdf("Boxplot_NegDeltaCqTrimmed_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(t(GabaGlu_NegDeltaCqTrimmed_AllSubjects_QCed)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCqTrimmed per sample for all genes", xlab="Sample ID", ylab="-DeltaCqTrimmed")
dev.off()


sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_QCed))
#[1] 133
dim(GabaGlu_NegDeltaCq_AllSubjects_QCed)
#[1] 142  96

GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed<-t(scale(GabaGlu_NegDeltaCq_AllSubjects_QCed, center=T, scale=T))
head(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)
dim(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)
#[1]  96 142

GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed<-t(scale(GabaGlu_NegDeltaCqTrimmed_AllSubjects_QCed, center=T, scale=T))
head(GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed)
dim(GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed)
#[1]  96 142

#Sanity check:
apply(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, 1, mean)
apply(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, 1, sd)
#Looks good.
apply(GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed, 1, mean)
apply(GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed, 1, sd)


sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed))
#[1] 133
#Since this is scale data, the mean value is 0 for all genes. Let's replace NAs with that:
GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed[is.na(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)]<-0

sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed))
#[1] 0

GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed[is.na(GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed)]<-0
sum(is.na(GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed))
#[1] 0


pdf("Boxplot_NegDeltaCqZscore_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of NegDeltaCq Z-scores per sample for all genes", xlab="Sample ID", ylab="NegDeltaCq Z-score")
dev.off()


#For comparison with earlier Cq plot:
pdf("Boxplot_NegDeltaCqZscore_perSample_AllGenes_ylim10.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of NegDeltaCq Z-scores per sample for all genes", xlab="Sample ID", ylab="NegDeltaCq Z-score", ylim=c(-10,10))
dev.off()


pdf("Boxplot_NegDeltaCqTrimmedZscore_perSample_AllGenes_ylim10.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of NegDeltaCqTrimmed Z-scores per sample for all genes", xlab="Sample ID", ylab="NegDeltaCqTrimmed Z-score", ylim=c(-10,10))
dev.off()
#I think the regular mean housekeeping gene correction actually performed slightly better.


#Let's look at the "housekeeping genes" specifically:

pdf("Boxplot_NegDeltaCqZscore_perSample_HousekeepingGenes_ylim10.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed[c(11, 86:96),]), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of NegDeltaCq Z-scores per sample for Housekeeping genes", xlab="Sample ID", ylab="NegDeltaCq Z-score", ylim=c(-10,10))
dev.off()

pdf("Boxplot_NegDeltaCqTrimmedZscore_perSample_HousekeepingGenes_ylim10.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCqTrimmed_AllSubjects_Scaled_QCed[c(11, 86:96),]), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of NegDeltaCqTrimmed Z-scores per sample for Housekeeping genes", xlab="Sample ID", ylab="NegDeltaCqTrimmed Z-score", ylim=c(-10,10))
dev.off()
#Ditto: I think the regular mean housekeeping gene correction actually performed slightly better.


#For comparison
pdf("Boxplot_CqZscore_perSample_HousekeepingGenes_ylim10.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_Scaled_QCed)[c(11, 86:96),], cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for Housekeeping genes", xlab="Sample ID", ylab="Cq Z-score", ylim=c(-10,10))
dev.off()
#Yes, the housekeeping genes are tightly variable across samples. They definitely seem to vary with card, but also vary by sample in a manner that does not necessarily correlate across replicates.


GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed<-cor(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, use="pairwise.complete.obs")
GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed

pdf("GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed.pdf", height=14, width=14)
heatmap(GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed)
dev.off()

#Trying PCA again:
pca_GabaGlu<-prcomp(t(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed))
tmp<-pca_GabaGlu$x[,1:10]
dim(tmp)
#[1] 142  10
rownames(tmp)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)
write.csv(tmp, "PCA_GabaGlu_NegDeltaCq.csv")

tmp<-pca_GabaGlu$rotation[,1:10]
write.csv(tmp, "pca_GabaGlu_NegDeltaCq_Eigenvectors.csv")

png("PCA_ScreePlot_GabaGlu_NegDeltaCq.png")
plot(summary(pca_GabaGlu)$importance[2,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 is still a huge amount of the variance (>30%) but much less than before (it was >60% before)

png("PCA_ScreePlot2_GabaGlu_NegDeltaCq.png")
plot(summary(pca_GabaGlu)$importance[3,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()


png("PC1vsPC2_GabaGlu_NegDeltaCq_byDiagnosis.png")
plot(pca_GabaGlu$x[,1]~pca_GabaGlu$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2", col=as.factor(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis))
dev.off()
#Not diagnosis - and there are 4 outliers in PC1

#Outliers:
pca_GabaGlu$x[,1][pca_GabaGlu$x[,1]>15]
#      37       49       13       16 
#28.62640 35.55227 19.57663 18.01706 

pca_GabaGlu$x[,1][pca_GabaGlu$x[,1]>20]
# 37       49 
# 28.62640 35.55227 

GabaGlu_OutlierSubjects<-pca_GabaGlu$x[,1]>20
sum(GabaGlu_OutlierSubjects)
#[1] 2

png("PC1_GabaGlu_NegDeltaCq_byCohort.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Cohort, xlab="Cohort", ylab="PC1")
dev.off()
#Some variation by cohort

pdf("PC1_GabaGlu_NegDeltaCq_byCard.pdf", width=10, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card, xlab="Cohort", ylab="PC1", cex.lab=0.3, las=3)
dev.off()
#Big variation by card.

summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card))
#esp cards 12, 19**, 25**,7
# Residual standard error: 5.034 on 106 degrees of freedom
# Multiple R-squared:  0.4636,	Adjusted R-squared:  0.2865 
# F-statistic: 2.617 on 35 and 106 DF,  p-value: 8.456e-05

pdf("PC1_GabaGlu_NegDeltaCq_byID.pdf", width=10, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID, ylim=c(-15, 45))
dev.off()


pdf("PC1_GabaGlu_NegDeltaCq_byID_CardEffectRemoved.pdf", width=10, height=5)
boxplot(summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card))$residuals~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID, ylim=c(-15, 45))
dev.off()
#Yeah, removing the card effect really tighten things up between replicate samples.

#Outputting a version where I can actually read the sample labels:
pdf("PC1_GabaGlu_NegDeltaCq_byID_CardEffectRemoved_longer.pdf", width=20, height=5)
boxplot(summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Card))$residuals~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID, ylim=c(-15, 45), cex.axis=0.5)
dev.off()
#Samples 37 and 49 look like the may still have some sort of bad measurement (very wide spread in their PC1 scores, one measurement has very high PC1)
#Note - these samples also have notably weird -Delta Cq distributions for all genes - the strangest in the experiment.
#As well as low houskeeping gene Cq z-scores.
#Samples 13 and 16 have high PC1, but their replicate measurements are actually reasonably consistent once card is controlled for.
#Weirdly, sample 1 looks consistent between replicates... despite having weird measurements in general.


pdf("PC1_GabaGlu_NegDeltaCq_byTrimmedMeanHousekeepingZScore.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~GabaGlu_TrimmedMeanHousekeepingZScore)
dev.off()
#The four outliers in PC1 have particularly low trimmed mean housekeeping gene expression.
#Otherwise, no correlation.


png("PC1_GabaGlu_NegDeltaCq_byDissection.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Dissecton.Group, xlab="Dissecton.Group", ylab="PC1")
dev.off()
#Not dissection

png("PC1_GabaGlu_NegDeltaCq_byGender.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Gender, xlab="Gender", ylab="PC1")
dev.off()
#Not gender

png("PC1_GabaGlu_NegDeltaCq_byAge.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Age, xlab="Age", ylab="PC1")
dev.off()
#Maybe a mild correlation with age - perhaps once the outliers are gone?

#Note: Subject #4463 is the one missing pH info - currently just has average pH

png("PC1_GabaGlu_byPH_NegDeltaCq.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$pH, xlab="pH", ylab="PC1")
dev.off()
#Maybe a mild correlation with pH - perhaps once the outliers are gone?

png("PC1_GabaGlu_NegDeltaCq_byHoursFinal.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Hours.Final, xlab="Hours.Final", ylab="PC1")
dev.off()
#Maybe a mild correlation with PMI - perhaps once the outliers are gone?

png("PC1_GabaGlu_NegDeltaCq_byID.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID, xlab="ID", ylab="PC1")
dev.off()
#Not ID number

png("PC1_GabaGlu_NegDeltaCq_bySubjectNumber.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Subject.Number, xlab="Subject.Number", ylab="PC1")
dev.off()
#Not subject number




#One more time: Removing the outlier samples:
pca_GabaGlu$x[,1][pca_GabaGlu$x[,1]>20]

dim(GabaGlu_NegDeltaCq_AllSubjects_QCed)
#[1] 142  96

GabaGlu_NegDeltaCq_AllSubjects_QCed2<-GabaGlu_NegDeltaCq_AllSubjects_QCed[pca_GabaGlu$x[,1]<20,]
dim(GabaGlu_NegDeltaCq_AllSubjects_QCed2)
#[1] 140  96

dim(SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
#[1] 142  17
SubjectInfo_OrderedForGabaGluCqMatrix_QCed2<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed[pca_GabaGlu$x[,1]<20,]
dim(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
#[1] 140  17

GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2<-t(scale(GabaGlu_NegDeltaCq_AllSubjects_QCed2, center=T, scale=T))
head(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2)
dim(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2)
#[1]  96 140

#Sanity check:
apply(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2, 1, mean)
apply(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2, 1, sd)


sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2))
#[1] 131
#Since this is scale data, the mean value is 0 for all genes. Let's replace NAs with that:
GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2[is.na(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2)]<-0

sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2))
#[1] 0



#Trying PCA again:
pca_GabaGlu<-prcomp(t(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2))
tmp<-pca_GabaGlu$x[,1:10]
dim(tmp)
#[1] 140  10
rownames(tmp)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed2)
write.csv(tmp, "PCA_GabaGlu_NegDeltaCq_NoOutliers.csv")

tmp<-pca_GabaGlu$rotation[,1:10]
write.csv(tmp, "pca_GabaGlu_NegDeltaCq_NoOutliers_Eigenvectors.csv")

png("PCA_ScreePlot_GabaGlu_NegDeltaCq_NoOutliers.png")
plot(summary(pca_GabaGlu)$importance[2,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 is a little below >30% of the variance still.

png("PCA_ScreePlot2_GabaGlu_NegDeltaCq_NoOutliers.png")
plot(summary(pca_GabaGlu)$importance[3,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()


png("PC1vsPC2_GabaGlu_NegDeltaCq_NoOutliers_byDiagnosis.png")
plot(pca_GabaGlu$x[,1]~pca_GabaGlu$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2", col=as.factor(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Diagnosis))
dev.off()
#Not diagnosis

png("PC1_GabaGlu_NegDeltaCq_NoOutliers_byCohort.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Cohort, xlab="Cohort", ylab="PC1")
dev.off()
#Some variation by cohort

pdf("PC1_GabaGlu_NegDeltaCq_NoOutliers_byCard.pdf", width=10, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, xlab="Cohort", ylab="PC1", cex.lab=0.3, las=3)
dev.off()
#Big variation by card.

summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card))
#esp cards 12, 22, 7  - not 19** and 25** anymore - they still have elevated PC1, but their sample size is lower.
# Residual standard error: 4.474 on 104 degrees of freedom
# Multiple R-squared:  0.4658,	Adjusted R-squared:  0.2861 
# F-statistic: 2.591 on 35 and 104 DF,  p-value: 0.0001051

pdf("PC1_GabaGlu_NegDeltaCq_NoOutliers_byID.pdf", width=10, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID, ylim=c(-15, 45))
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_NoOutliers_byID_CardEffectRemoved.pdf", width=10, height=5)
boxplot(summary.lm(lm(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card))$residuals~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID, ylim=c(-15, 45))
dev.off()
#Yeah, removing the card effect really tighten things up between replicate samples.

pdf("PC1_GabaGlu_NegDeltaCq_NoOutliers_byTrimmedMeanHousekeepingZScore.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~GabaGlu_TrimmedMeanHousekeepingZScore[OutlierSubjects==F])
dev.off()
# #There is an outlier again with a PC1 that is high and particularly high housekeeping gene expression.
GabaGlu_TrimmedMeanHousekeepingZScore[OutlierSubjects==F]>3
# #It is the second sample 12
#Not including that subject, there semes to be a negative correlation with housekeeping z-score:
cor(pca_GabaGlu$x[,1],GabaGlu_TrimmedMeanHousekeepingZScore[OutlierSubjects==F])
#[1] -0.3237057

png("PC1_GabaGlu_NegDeltaCq_NoOutliers_byDissection.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group, xlab="Dissection.Group", ylab="PC1")
dev.off()
#Not dissection

png("PC1_GabaGlu_NegDeltaCq_NoOutliers_byGender.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Gender, xlab="Gender", ylab="PC1")
dev.off()
#Not gender

png("PC1_GabaGlu_NegDeltaCq_NoOutliers_byAge.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Age, xlab="Age", ylab="PC1")
dev.off()
#Maybe a mild correlation with age 
cor(pca_GabaGlu$x[,1],SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Age)
#[1] 0.2758202

#Note: Subject #4463 is the one missing pH info - currently just has average pH

png("PC1_GabaGlu_byPH_NoOutliers_NegDeltaCq.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$pH, xlab="pH", ylab="PC1")
dev.off()
#Maybe a mild negative correlation with pH

png("PC1_GabaGlu_NegDeltaCq_NoOutliers_byHoursFinal.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Hours.Final, xlab="Hours.Final", ylab="PC1")
dev.off()
#Maybe a mild correlation with PMI 


png("PC1_GabaGlu_NegDeltaCq_NoOutliers_byID.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID, xlab="ID", ylab="PC1")
dev.off()
#Not ID number

png("PC1_GabaGlu_NegDeltaCq_NoOutliers_bySubjectNumber.png")
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Subject.Number, xlab="Subject.Number", ylab="PC1")
dev.off()
#Not subject number but maybe a cohort effect.


#Formally examining the relationship between PC1 and the usual suspects & batch effects:

Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)

#With card as a fixed effect (so that we can easily see how much it matters)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Cohort+Card + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)       3.6205  1  0.0570727 .  
# Diagnosis         9.1333  2  0.0103929 *  
#   Age              26.3019  1   2.92e-07 ***
#   Gender            1.6052  1  0.2051666    
# pH                5.2608  1  0.0218109 *  
#   Hours.Final      11.8265  1  0.0005839 ***
#   Dissecton.Group  16.8457 11  0.1125241    
# Cohort           23.3694  8  0.0029210 ** 
#   Card            283.5835 31  < 2.2e-16 ***

#With card as a random effect (probably more appropriate)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Cohort+ (1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)      2.7538  1   0.097025 .  
# Diagnosis        5.9543  2   0.050938 .  
# Age             22.1046  1  2.582e-06 ***
#   Gender           0.2790  1   0.597332    
# pH               4.1956  1   0.040528 *  
#   Hours.Final      9.9671  1   0.001594 ** 
#   Dissecton.Group  3.3443 11   0.985349    
# Cohort          13.5930  8   0.093012 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Looks like Card, Age, pH, PMI and cohort really matter as co-variates (no surprise)
#Most of the variables seem to have stronger effects when Card is modeled as a fixed effect (e.g., Diagnosis). I wonder why?
#Dissection group doesn't seem to matter as much - maybe we can toss it out as a co-variate?

#With card and dissection group as a random effect (probably more appropriate)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+ (1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)  2.6875  1   0.101141    
# Diagnosis    5.3499  2   0.068909 .  
# Age         20.7949  1  5.112e-06 ***
#   Gender       0.6418  1   0.423045    
# pH           3.8154  1   0.050783 .  
# Hours.Final  9.5396  1   0.002011 ** 
#   Cohort      12.3364  8   0.136815    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting - including dissecton group as a random effect seems to weaken all of the other effects even more. 


#Formally examining the relationship between PC2 and the usual suspects & batch effects:

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)

#With card as a fixed effect (so that we can easily see how much it matters)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Cohort+Card + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)       5.0396  1    0.02477 *  
#   Diagnosis         3.2239  2    0.19949    
# Age               2.8012  1    0.09419 .  
# Gender            4.2924  1    0.03828 *  
#   pH                4.1574  1    0.04145 *  
#   Hours.Final       0.9137  1    0.33913    
# Dissecton.Group  12.7244 11    0.31172    
# Cohort           23.8910  8    0.00239 ** 
#   Card            111.1129 31   5.94e-11 ***

#With card as a random effect
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Cohort+ (1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept)      4.1020  1    0.04283 *
#   Diagnosis        0.4175  2    0.81159  
# Age              0.2676  1    0.60497  
# Gender           0.2069  1    0.64917  
# pH               3.5991  1    0.05781 .
# Hours.Final      2.5801  1    0.10822  
# Dissecton.Group 22.5137 11    0.02068 *
#   Cohort          16.9747  8    0.03037 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Card really matters (again)
#Dissection group kind-of matters for PC2. So does cohort, pH and gender.

#Hmm.... I guess we should include all of them.

#With card and dissecton group as a random effect
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+ (1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)   
# (Intercept)  6.6856  1    0.00972 **
#   Diagnosis    0.3585  2    0.83589   
# Age          0.3349  1    0.56277   
# Gender       0.2867  1    0.59233   
# pH           6.4296  1    0.01122 * 
#   Hours.Final  1.3591  1    0.24369   
# Cohort      13.3646  8    0.09991 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Let's see if that reduced the batch effects:
#Housekeeping genes first (even though it is a little circular):
for(i in c(11, 86:96)){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Cohort+ (1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndBatchFixed_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndBatchFixed_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#Nope - the effects of cohort and dissection group are still there, and now some of the other effects are stronger (Gender, pH). Two housekeeping genes have diagnosis effects: HPRT1 and TFRC (HMBS is a trend).


#How does this compare with the results from the average housekeeping gene signal?
#With the outliers:
Temp<-data.frame(y=(GabaGlu_TrimmedMeanHousekeeping*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+Card + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 97.7940  1  < 2.2e-16 ***
#   Diagnosis    2.2044  2   0.332136    
# Age          2.3057  1   0.128900    
# Gender       0.1408  1   0.707448    
# pH           0.2082  1   0.648204    
# Hours.Final  4.3401  1   0.037226 *  
#   Cohort      23.4563  8   0.002825 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Without the outliers:
Temp<-data.frame(y=(GabaGlu_TrimmedMeanHousekeeping[(OutlierSubjects==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 115.7952  1  < 2.2e-16 ***
#   Diagnosis     3.4498  2   0.178189    
# Age           2.5146  1   0.112797    
# Gender        0.0143  1   0.904698    
# pH            0.0074  1   0.931321    
# Hours.Final   5.3417  1   0.020821 *  
#   Cohort       20.7965  8   0.007708 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Housekeeping genes vs. just diagnosis, no outliers
Model<-lmer(y~Diagnosis + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 1.0229e+05  1     <2e-16 ***
#   Diagnosis   1.5537e+00  2     0.4598    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Vs. just diagnosis and biological co-variates (no batch):
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 158.2512  1    < 2e-16 ***
#   Diagnosis     2.5160  2    0.28422    
# Age           3.0992  1    0.07833 .  
# Gender        0.0011  1    0.97333    
# pH            0.9918  1    0.31930    
# Hours.Final   3.4225  1    0.06431 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Vs. just diagnosis, biological co-variates, and cohort (no dissection batch or card):
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 132.8228  1  < 2.2e-16 ***
#   Diagnosis     4.2890  2   0.117127    
# Age           1.0379  1   0.308315    
# Gender        0.1185  1   0.730676    
# pH            0.4862  1   0.485620    
# Hours.Final   6.7073  1   0.009602 ** 
#   Cohort       15.0218  8   0.058722 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     130.9883  1    < 2e-16 ***
#   Diagnosis         3.7184  2    0.15580    
# Age               1.3664  1    0.24243    
# Gender            0.0135  1    0.90735    
# pH                0.0731  1    0.78686    
# Hours.Final       3.2327  1    0.07218 .  
# Cohort           16.2303  8    0.03920 *  
#   Dissecton.Group  21.5115 11    0.02844 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

#fixed-effect model matrix is rank deficient so dropping 4 columns / coefficients
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     118.4116  1  < 2.2e-16 ***
#   Diagnosis         1.3611  2  0.5063430    
# Age              10.5804  1  0.0011430 ** 
#   Gender            1.8888  1  0.1693358    
# pH                1.5465  1  0.2136563    
# Hours.Final       6.4484  1  0.0111050 *  
#   Cohort           26.2469  8  0.0009528 ***
#   Dissecton.Group  29.9598 11  0.0016078 ** 
#   Card            161.5681 31  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Looking at regular mean housekeeping (not trimmed mean):
Temp<-data.frame(y=(GabaGlu_MeanHousekeeping[(OutlierSubjects==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 103.7362  1    < 2e-16 ***
#   Diagnosis     2.8674  2    0.23842    
# Age           2.0019  1    0.15710    
# Gender        0.0017  1    0.96749    
# pH            0.0008  1    0.97773    
# Hours.Final   4.3467  1    0.03708 *  
#   Cohort       20.0211  8    0.01026 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Housekeeping genes vs. just diagnosis, no outliers
Model<-lmer(y~Diagnosis + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 93591.0104  1     <2e-16 ***
#   Diagnosis       1.1219  2     0.5707    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Vs. just diagnosis and biological co-variates (no batch):
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 150.2023  1    < 2e-16 ***
#   Diagnosis     2.0675  2    0.35568    
# Age           2.2393  1    0.13454    
# Gender        0.0008  1    0.97748    
# pH            1.5246  1    0.21692    
# Hours.Final   3.0098  1    0.08276 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 116.4929  1     <2e-16 ***
#   Diagnosis     3.5545  2     0.1691    
# Age           0.5751  1     0.4483    
# Gender        0.2984  1     0.5849    
# pH            0.8825  1     0.3475    
# Hours.Final   5.6758  1     0.0172 *  
#   Cohort       12.0540  8     0.1488    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     110.3844  1    < 2e-16 ***
#   Diagnosis         2.7959  2    0.24711    
# Age               0.7350  1    0.39128    
# Gender            0.0370  1    0.84743    
# pH                0.3674  1    0.54444    
# Hours.Final       2.0358  1    0.15363    
# Cohort           12.1745  8    0.14359    
# Dissecton.Group  19.1991 11    0.05761 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

#fixed-effect model matrix is rank deficient so dropping 4 columns / coefficients
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     110.7631  1  < 2.2e-16 ***
#   Diagnosis         1.1582  2  0.5604035    
# Age               7.4788  1  0.0062429 ** 
#   Gender            1.3980  1  0.2370505    
# pH                0.9870  1  0.3204655    
# Hours.Final       5.6705  1  0.0172526 *  
#   Cohort           26.3371  8  0.0009195 ***
#   Dissecton.Group  34.0645 11  0.0003530 ***
#   Card            241.2935 31  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Which cohorts, dissection groups, and cards are causing the problems?
summary(lmer(y~Age+Hours.Final+Cohort+Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F))
#Dissection group6, dissection group 8, dissection group 11 - maybe we could cluster these instead?
#Card 10, 17, 2, 22, 31-34, 36, 6, 7
#Cohort 12, Cohort 6, Cohort 2

library(dplyr)

rescov <- function(model, data) {
  var.d <- crossprod(getME(model,"Lambdat"))
  Zt <- getME(model,"Zt")
  vr <- sigma(model)^2
  var.b <- vr*(t(Zt) %*% var.d %*% Zt)
  sI <- vr * Diagonal(nrow(data))
  var.y <- var.b + sI
  invisible(var.y)
}


mod1<-lmer(y~Cohort+(1 | ID), data = Temp  %>% arrange(Cohort), REML=F)
mod1<-lmer(y~1+(1 | ID), data = Temp  %>% arrange(Cohort), REML=F)
#mod1 <- lmer(diameter ~ 1 + (1 | sample), data = Penicillin %>% arrange(sample))
rc1 <- rescov(mod1, Temp)
image(rc1)
#Well that was unhelpful.

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+ (1 | ID), data = Temp %>% arrange(Cohort), REML=F)
plot(summary(Model)$residuals)
boxplot(summary(Model)$residuals~(Temp %>% arrange(Cohort))$Cohort)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+ (1 | ID), data = Temp %>% arrange(Card), REML=F)
plot(summary(Model)$residuals)
boxplot(summary(Model)$residuals~(Temp %>% arrange(Card))$Card)



#Simplest model: Just the Diagnosis effects (without controlling for the biggest co-variates)
for(i in c(1:96)){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+ (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByJustDiagnosis_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByJustDiagnosis_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}


# Model with biological co-variates but not batch variables:
for(i in c(1:96)){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+ (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspects_NoBatch_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspects_NoBatch_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}




for(i in c(1:96)){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}



#How does this compare with the results if we don't subtract out average housekeeping gene signal?

for(i in c(1:96)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_Cq.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_Cq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}


#I'm a little worried about the lack of effect in SST & PVALB. I'm worried that the cohorts are collinear with diagnosis and somehow subtracting out the effect.
#Let's try a version without cohort.
for(i in c(1:96)){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq_NoCohort.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq_NoCohort.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#SST has a diagnosis effect if Cohort isn't included in the model.
#But many of the age effects weakened and some of the diagnosis effects (although some were also gained). Actually most everything else seems weaker, which makes sense since PC1 related to cohort.
#And 18S relates to diagnosis then - which is a little worrisome.

#Dissection group is not related to PC1 and uses up a whole bunch of degrees of freedom. What if we yank it instead?

for(i in c(1:96)){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq_NoDissection.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq_NoDissection.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#That output is identical to the output with dissection.... why??? 


#Let's try comparing the version with dissection as a fixed effect
for(i in c(1:96)){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+(1 | Card) + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq_Dissection.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq_Dissection.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Dissection is significant for quite a few of the genes...I guess it might be hard to cut it out.


#Is modeling things via random effect just causing issues of some sort?
for(i in c(1:96)){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+Card + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq_DissectionCard.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndBatchRandom_NegDeltaCq_DissectionCard.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#"fixed-effect model matrix is rank deficient so dropping 4 columns / coefficients"
#Makes a big difference... but is it legitimate??? Definitely getting close to the number of variables (57 Df) as there are subjects :(
#but I'm definitely worried about how card is being handled as a random effect variable. 


table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID)

# 1 10 12 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31 32 33 34 37 38 39  4 40 41 42 43 44 45 46 47 48 49 
# 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  1  2  2  2  2  2  2  2  2  1  1  2  2  2  2  2  2  2  2  2  2  2  1 
# 5 50 51 52 53 54 55 56 57 58 59  6 60 61 62 63 64 65 66 67 68 69  7 70 71 72 73 74 75  8  9 
# 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 

dim(GabaGlu_NegDeltaCq_AllSubjects_QCed2)
#[1] 140  96
head(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

GabaGlu_NegDeltaCq_AllSubjects_QCed2_AveByID<-matrix(0, 72, 96)
row.names(GabaGlu_NegDeltaCq_AllSubjects_QCed2_AveByID)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID))
colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2_AveByID)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)
head(GabaGlu_NegDeltaCq_AllSubjects_QCed2_AveByID)

for(i in c(1:96)){
GabaGlu_NegDeltaCq_AllSubjects_QCed2_AveByID[,i]<-tapply(GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID, function(y) mean(y, na.rm=T))
}

head(GabaGlu_NegDeltaCq_AllSubjects_QCed2_AveByID)


#Adriana was worried about diagnosis effects or noise in the housekeeping genes interfering with the ability to detect diagnosis effects - what if we controlled for housekeeping signal in a larger model that included diagnosis istead of just normalizing by it?

GabaGlu_Cq_AllSubjects_QCed[(OutlierSubjects==F)] GabaGlu_MeanHousekeeping[(OutlierSubjects==F)]

Temp<-data.frame(y=(GabaGlu_Cq_AllSubjects_QCed[i, (OutlierSubjects==F)]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)

# Temp<-data.frame(y=(GabaGlu_Cq_AllSubjects_QCed[i, (OutlierSubjects==F)]*-1), HK=(GabaGlu_TrimmedMeanHousekeeping[(OutlierSubjects==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)

plot(Temp$y~Temp$HK)

#Which subject is it that has the crazy low HK signal?

Temp[Temp$HK<(-23),]

# y        HK ID       Card Subject.Number  Barcode    Cohort Diagnosis Age Gender   pH AFS Hours.Final
# 12577 -22.196 -24.03392 12 Card 6.eds           4830 B007529A Cohort 11   Control  53      M 6.77   0       23.25
# Slab.Format Slab.Number Dissecton.Group SampleNumber QC_AmplificationIssue CohortDissection
# 12577        SLAB           1               2            3                     N      Cohort 11.2

#Yeah, going back to the Cq Z scores, that sample just has really high Cq for everything, and once it is corrected (delta Cq) doesn't seem particularly odd - at least, in terms of PC1.
#It's going to really skew any results using HK as co-variate though - huge leverage.

OutlierSubjects2<-OutlierSubjects

OutlierSubjects2[GabaGlu_MeanHousekeeping>(24)]
# 12 
# FALSE 
OutlierSubjects2[GabaGlu_MeanHousekeeping>(24)]<-TRUE

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed[(OutlierSubjects2==F),]
dim(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
#[1] 139  17

# Temp<-data.frame(y=(GabaGlu_Cq_AllSubjects_QCed[i, (OutlierSubjects2==F)]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
# 
# Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+HK+ (1 | ID), data = Temp, REML=F)
# car::Anova(Model, type="III")

#Does that change the relationship between housekeeping genes and other variables at all?
#Looking at regular mean housekeeping (not trimmed mean):
Temp<-data.frame(y=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+(1|Dissecton.Group)+(1 | Card) + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 178.3475  1  < 2.2e-16 ***
#   Diagnosis     3.9335  2   0.139908    
# Age           2.3424  1   0.125899    
# Gender        0.0656  1   0.797898    
# pH            0.0005  1   0.981947    
# Hours.Final   5.6491  1   0.017464 *  
#   Cohort       26.8717  8   0.000744 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Cohort, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis)
#Cohort 6, 7 & 8 are all BP and Schiz. :(  

#                 Control BP Schiz
# Cohort 11           11  2     6
# Cohort 12           17  8     4
# Cohort 13            9  0     4
# Cohort 7             0 11    11
# Cohort 8             0  4     8
# Dep Cohort 1         4  4     0
# Dep Cohort 5         2  8     0
# Dep Cohort 6         0  4     4
# Schiz Cohort 2       8  0    10

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Cohort, SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis)
#That was also true before removing the outlier subjects


#Housekeeping genes vs. just diagnosis, no outliers
Model<-lmer(y~Diagnosis + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 123889.906  1     <2e-16 ***
#   Diagnosis        2.367  2     0.3062    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Vs. just diagnosis and biological co-variates (no batch):
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 190.0670  1    < 2e-16 ***
#   Diagnosis     3.2552  2    0.19640    
# Age           3.1051  1    0.07805 .  
# Gender        0.0122  1    0.91219    
# pH            1.3768  1    0.24064    
# Hours.Final   3.8402  1    0.05004 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 142.5448  1    < 2e-16 ***
#   Diagnosis     4.0248  2    0.13367    
# Age           0.9247  1    0.33625    
# Gender        0.3602  1    0.54838    
# pH            0.6503  1    0.42000    
# Hours.Final   6.0339  1    0.01403 *  
#   Cohort       12.0920  8    0.14715    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     145.1362  1    < 2e-16 ***
#   Diagnosis         3.4864  2    0.17496    
# Age               0.7052  1    0.40103    
# Gender            0.2578  1    0.61166    
# pH                0.5135  1    0.47361    
# Hours.Final       2.5891  1    0.10760    
# Cohort           16.4577  8    0.03628 *  
#   Dissecton.Group  22.2724 11    0.02235 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Cohort+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     197.7704  1    < 2e-16 ***
#   Diagnosis         2.7593  2    0.25167    
# Age               3.5817  1    0.05842 .  
# Gender            0.0426  1    0.83653    
# pH                1.5323  1    0.21576    
# Hours.Final       1.9175  1    0.16614    
# Dissecton.Group  17.7434 11    0.08773 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
#fixed-effect model matrix is rank deficient so dropping 4 columns / coefficients

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     198.9689  1  < 2.2e-16 ***
#   Diagnosis         4.7442  2   0.093284 .  
# Age               7.0715  1   0.007832 ** 
#   Gender            1.1520  1   0.283126    
# pH                0.0059  1   0.938819    
# Hours.Final       0.9707  1   0.324506    
# Dissecton.Group  50.0656 11  6.092e-07 ***
#   Card            252.6398 31  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Card+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 172.6900  1    < 2e-16 ***
#   Diagnosis     1.8080  2    0.40494    
# Age           6.3044  1    0.01204 *  
#   Gender        0.2447  1    0.62085    
# pH            0.2725  1    0.60168    
# Hours.Final   0.8061  1    0.36926    
# Card        251.4338 35    < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



Temp<-data.frame(y=(GabaGlu_Cq_AllSubjects_QCed[i, (OutlierSubjects2==F)]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

for(i in c(1:96)){
  Temp<-data.frame(y=(GabaGlu_Cq_AllSubjects_QCed[i, (OutlierSubjects2==F)]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+HK+(1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#That is so much prettier now. :)

#For direct comparison:
dim(GabaGlu_NegDeltaCq_AllSubjects_QCed)


for(i in c(1:96)){
  Temp<-data.frame(y=(GabaGlu_NegDeltaCq_AllSubjects_QCed[(OutlierSubjects2==F), i]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+(1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Very similar results, but some of the relationships are a little stronger when using HK as a co-variate.



for(i in c(1:96)){
  Temp<-data.frame(y=(GabaGlu_Cq_AllSubjects_QCed[i, (OutlierSubjects2==F)]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+HK+(1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByJustDiagnosis_NegCq_HKcov.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByJustDiagnosis__NegCq_HKcov.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#For the most part you can't see the effects unless you add a few co-variates to clean it up.


#Let's output things better so we can add FDR, etc, plot trends, etc.

Temp<-data.frame(y=(GabaGlu_Cq_AllSubjects_QCed[i, (OutlierSubjects2==F)]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+HK+(1 | ID), data = Temp, REML=F)
names(car::Anova(Model, type="III"))

ChiSquare<-car::Anova(Model, type="III")[[1]]

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare<-matrix(0, 96, 8)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval<-matrix(0, 96, 8)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)


for(i in c(1:96)){
  Temp<-data.frame(y=(GabaGlu_Cq_AllSubjects_QCed[i, (OutlierSubjects2==F)]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+HK+(1 | ID), data = Temp, REML=F)
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}



Temp<-data.frame(y=(GabaGlu_NegDeltaCq_AllSubjects_QCed[(OutlierSubjects2==F), i]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+(1 | ID), data = Temp, REML=F)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare<-matrix(0, 96, 7)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval<-matrix(0, 96, 7)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed)


for(i in c(1:96)){
  Temp<-data.frame(y=(GabaGlu_NegDeltaCq_AllSubjects_QCed[(OutlierSubjects2==F), i]*-1), HK=(GabaGlu_MeanHousekeeping[(OutlierSubjects2==F)]*-1), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+(1 | ID), data = Temp, REML=F)
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

#Alright - now we can directly compare the two methods:

colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare)
plot(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare[,2]~MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare[,2] )
summary.lm(lm(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare[,2]~MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare[,2]))
# Call:
#   lm(formula = MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare[, 
#                                                                                          2] ~ MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare[, 
#                                                                                                                                                                        2])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.3164 -0.5251 -0.1119  0.4649  2.8981 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                    0.59238    0.13298   4.455 2.31e-05 ***
#   MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare[, 2]  0.95815    0.03775  25.380  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9446 on 94 degrees of freedom
# Multiple R-squared:  0.8727,	Adjusted R-squared:  0.8713 
# F-statistic: 644.2 on 1 and 94 DF,  p-value: < 2.2e-16

#Surprisingly, the Chi-Square values for diagnosis are actually a little smaller in the HKcov group.
#But in general the results are strongly correlated.


plot(log10(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval[,2])~log10(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[,2]))
summary.lm(lm(log10(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval[,2])~log10(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[,2])))

# Call:
#   lm(formula = log10(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval[, 
#                                                                                           2]) ~ log10(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[, 
#                                                                                                                                                                           2]))
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.62931 -0.10095  0.02431  0.11402  0.72014 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                     -0.12863    0.02888  -4.455 2.31e-05 ***
#   log10(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[, 2])  0.95815    0.03775  25.380  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Residual standard error: 0.2051 on 94 degrees of freedom
# Multiple R-squared:  0.8727,	Adjusted R-squared:  0.8713 
# F-statistic: 644.2 on 1 and 94 DF,  p-value: < 2.2e-16

#Either way, there are only a handful of nominally significant effects of diagnosis:


MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval[,2]<0.05,2]
# ABAT          SST         BDNF          18S       GABRR1        GRIK1       HOMER1         GFAP        MAPK1 
# 0.0200350952 0.0067581267 0.0425843277 0.0211890995 0.0437217580 0.0490498516 0.0060149755 0.0158131778 0.0185572312 
# SLC1A1        RPLP0         TFRC 
# 0.0325898381 0.0292360298 0.0005320054 

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[,2]<0.05,2]
# ABAT          SST       GABRR1       HOMER1         GFAP         TFRC 
# 0.0242785157 0.0067023957 0.0092038502 0.0297352912 0.0150597962 0.0002592586 

#Let's add FDR correction:


library(multtest)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval
MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval

for(i in c(1:8)){
tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval[,i], proc="BH")
MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)
}
head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR)

for(i in c(1:7)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[,i], proc="BH")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}
head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR[,2]<0.05,2]
#numeric(0)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR[,2]<0.10,2]
#[1] 0.05107252
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR[,2]<0.10]
#[1] "TFRC"

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR[,2]<0.05,2]
#[1] 0.02488883
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR[,2]<0.05]
#[1] "TFRC"

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR[,2]<0.10,2]
#[1] 0.02488883
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR[,2]<0.10]
#[1] "TFRC"

#Well... that's awkward. :(

#Speaking of which: I suppose it might be more appropriate to remove the housekeeping genes from the analysis before running multiple comparison corrections.


MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare_NoHK<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_ChiSquare[-c(11, 86:96),]
MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare_NoHK<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_ChiSquare[-c(11, 86:96),]

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval_NoHK<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval[-c(11, 86:96),]
MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval[-c(11, 86:96),]

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK
MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval_NoHK

for(i in c(1:8)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,i], proc="BH")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}
head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK)

for(i in c(1:7)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_Pval_NoHK[,i], proc="BH")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}
head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK)


MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK[,2]<0.10,2]
#numeric(0)
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK[,2]<0.10]
#character(0)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK[,2]<0.20,2]
#numeric(0)
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK[,2]<0.20]
#character(0)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK[,2]<0.10,2]
#numeric(0)
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK[,2]<0.10]
#character(0)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK[,2]<0.20,2]
#numeric(0)
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegDeltaCq_FDR_NoHK[,2]<0.20]
#character(0)

colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK)
hist(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,2], breaks=20)
hist(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,3], breaks=20)
#lots of Age-related genes
hist(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,4], breaks=20)
#Almost no gender-related genes
hist(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,5], breaks=20)
#lots of pH-related genes
hist(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,6], breaks=20)
#lots of PMI-related genes
hist(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,7], breaks=20)
#lots of dissection-related genes
hist(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,8], breaks=20)
#And basically everything is related to HK

#Is anything related to Gender? i.e., could we potentially throw it out?

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_FDR_NoHK[,4]<0.20]
#character(0)
#Not after correcting for FDR

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndDissection_NegCq_HKcov_Pval_NoHK[,4]<0.05,4]
# ADCY7      GABRD     GABRR2       GAD1       GRM4    SLC17A6    SLC6A11 
# 0.01048353 0.00958328 0.01821042 0.03796374 0.03628449 0.03372139 0.03835678 
#There are a few with nominal relationships

#To Do (since we need to organize and re-run everything anyway)
#Output mean housekeeping gene histogram
#Re-run pca after removing all 3 outliers
#Output diagnosis vs. all covariates *following outlier removal* to determine which are multicollinear.
#Maybe output a version of the analysis without gender included? (since it is somewhat multicollinear but potentially mostly meaningless? i.e. adding noise)

#Dealing with card:

#I'm hesitant to remove the effect of card completely because it also hypothetically includes 1/4 of the variation associated with each sample.
#Averaging the replicates would leave some of that card-related variation, but protect the sample-related variation.
#However, due to technical issues and outlier removal, not all subjects have replicate samples (5), and therefore are way more susceptible to card effects
#We could try doing some sort of modified averaging by sample - e.g., calculate the average across the two cards containing replicates for each gene, and then add that back to the residuals and average the replicates.

