#Example code for importing Gemma's differential expression results for our screened datasets
#Megan Hagenauer
#2026-07-22

#This is the code if we already narrowed our result sets down to only the ones that we want
ResultSet_contrasts<-read.csv("ResultSets_Screened.csv", header=TRUE, stringsAsFactors = FALSE )

#Here is code if we want to subset our result sets by the "Include" column:
list.files()
ResultSets_toScreen<-read.csv("ResultSets_toScreen - ResultSets_toScreen.csv", header=TRUE, stringsAsFactors = FALSE)

str(ResultSets_toScreen)

ResultSet_contrasts<-ResultSets_toScreen[ResultSets_toScreen$Include=="Y", ]

str(ResultSet_contrasts)

#Download the function from our github repository:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_DownloadingDEResults.R

#Reading in the function:

source("Function_DownloadingDEResults.R")

library(gemma.R)
library(tidyr)

DownloadingDEResults(ResultSet_contrasts)


#################


SavingGemmaDEResults_forEachResultSet<-function(differentials, UniqueResultSetIDs, ResultSet_contrasts){
  
  for (i in c(1:length(differentials))){
    
    ThisResultSet<-UniqueResultSetIDs[i]
    
    #Pulling out the dataset name from our other data frame
    #For some reason I can find this in the Gemma ResultSet differential expression output
    #Since some datasets have multiple result sets, we just grab the dataset name from the first entry
    ThisDataSet<-ResultSet_contrasts$ExperimentID[ResultSet_contrasts$ResultSetIDs==ThisResultSet][1] 
    
    #Write out a data frame containing the differential expression output for the result set
    #And name it with the dataset id and result set id:
    write.csv(differentials[[i]], paste("DEResults", ThisDataSet, ThisResultSet, ".csv", sep="_"))
    
    rm(ThisDataSet, ThisResultSet)
  }
}


#You can input this function by running the code discussed above to create the function in your R environment. 

#Alternatively, you can download the script for the function from our Github site and save the file in your working directory:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_SavingGemmaDEResults_forEachResultSet.R

#And then source it from your working directory:
source("Function_SavingGemmaDEResults_forEachResultSet.R")

#Example usage:

SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)
