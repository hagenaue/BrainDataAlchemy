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

#Reading in the function:

DownloadingDEResults<-function(ResultSet_contrasts){
  
  #Some ResultSets have more than one statistical contrast, so they are present more than once in our data frame, e.g.:
  #ResultSet_contrasts$ResultSetIDs
  #[1] 553805 553805 553805 570552 556647
  
  #To pull down the statistical results, we'll only want the unique result set ids:
  UniqueResultSetIDs<-unique(ResultSet_contrasts$ResultSetIDs)
  
  print("These are the result sets that you identified as being of interest")
  print(UniqueResultSetIDs)
  #553805 570552 556647
  
  differentials <- UniqueResultSetIDs %>% lapply(function(x){
    #   # take the first and only element of the output. the function returns a list 
    #   # because single experiments may have multiple resultSets. Here we use the 
    #   # resultSet argument to directly access the results we need
    get_differential_expression_values(resultSet = x)[[1]]
  })
  
  str(differentials)
  #That code worked. Excellent!
  
  # # some datasets might not have all the advertised differential expression results
  # # calculated due to a variety of factors. here we remove the empty differentials
  missing_contrasts <- differentials %>% sapply(nrow) %>% {.==0}
  #[1] FALSE FALSE FALSE
  differentials <<- differentials[!missing_contrasts]
  UniqueResultSetIDs<<-UniqueResultSetIDs[!missing_contrasts]
  
  print("These are the result sets that had differential expression results:")
  print(UniqueResultSetIDs)
  
  print("Your differential expression results for each of your result sets are stored in the object named differentials. This object is structured as a list of data frames. Each element in the list represetns a result set, with the data frame containing the differential expression results")
  
  #Within any particular Result Set, there are likely to be some contrasts that we want and others that we don't want
  #For example, a result set might contain a variety of stress interventions
  #And maybe we only want the acute stress contrast results
  
  #We already identified which statistical contrasts we wanted during our screening:
  #This is the object with the specific contrast ids that we want:
  #ResultSet_contrasts$ContrastIDs
  
  #Which will be these columns within the listed dataframes of differential expression results:
  
  print("These are the columns for the effect sizes for our statistical contrasts of interest (Log(2) Fold Changes")
  Contrasts_Log2FC<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_log2fc", sep="")
  
  print(Contrasts_Log2FC)
  #[1] "contrast_151618_log2fc" "contrast_151617_log2fc" "contrast_151619_log2fc"
  #[4] "contrast_186753_log2fc" "contrast_204289_log2fc"
  
  print("these are the columns for the T-statistics for our statistical contrasts of interest - we will use that information to derive the sampling variances")
  
  Contrasts_Tstat<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_tstat", sep="")
  
  print(Contrasts_Tstat)
  # [1] "contrast_151618_tstat" "contrast_151617_tstat" "contrast_151619_tstat"
  # [4] "contrast_186753_tstat" "contrast_204289_tstat"  
  
}


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
