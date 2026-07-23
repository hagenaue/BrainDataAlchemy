#Example code for processing the differential expression results for each study
#Megan Hagenauer
#2026-07-23

#####################

library(gemma.R)

install.packages("tidyverse")

library(tidyverse)

install.packages("dplyr")

library (dplyr)

########################

str(differentials)

str(differentials[[1]])

#Here is an example of pulling out the results of the first result set in the differentials object:

DE_Results<-differentials[[1]]

#The structure of the new DE_Results object:
str(DE_Results)
#Classes ‘data.table’ and 'data.frame':	17064 obs. of  12 variables:


sum(c(1,4))
sum(c(1,4,3))

sum(c(TRUE, FALSE, TRUE, FALSE))

# & AND
# | OR
# == equivalent
# != not equivalent

FilteringDEResults_GoodAnnotation<-function(DE_Results){
  
  print("# of rows in results")
  print(nrow(DE_Results))
  
  print("# of rows with missing NCBI annotation:")
  print(sum(DE_Results$NCBIid==""|DE_Results$NCBIid=="null"))
  
  print("# of rows with NA NCBI annotation:")
  print(sum(is.na(DE_Results$NCBIid)))
  
  print("# of rows with missing Gene Symbol annotation:")
  print(sum(DE_Results$GeneSymbol==""|DE_Results$GeneSymbol=="null"))
  
  print("# of rows mapped to multiple NCBI_IDs:")
  print(length(grep('\\|', DE_Results$NCBIid)))
  
  print("# of rows mapped to multiple Gene Symbols:")
  print(length(grep('\\|', DE_Results$GeneSymbol)))
  
  #I only want the subset of data which contains rows that do not contain an NCBI EntrezID of ""
  DE_Results_NoNA<-DE_Results[(DE_Results$NCBIid==""|DE_Results$NCBIid=="null")==FALSE & is.na(DE_Results$NCBIid)==FALSE,]
  
  #I also only want the subset of data that is annotated with a single gene (not ambiguously mapped to more than one gene)
  if(length(grep('\\|', DE_Results_NoNA$NCBIid))==0){
    DE_Results_GoodAnnotation<<-DE_Results_NoNA
  }else{
    #I only want rows annotated with a single Gene Symbol (no pipe):
    DE_Results_GoodAnnotation<<-DE_Results_NoNA[-(grep('\\|', DE_Results_NoNA$NCBIid)),]
  }
  #I used a double arrow in that conditional to place DE_Results_GoodAnnotation back out in the environment outside the function 
  
  print("# of rows with good annotation")
  print(nrow(DE_Results_GoodAnnotation))
  
  #For record keeping (sometimes useful for troubleshooting later)
  write.csv(DE_Results_GoodAnnotation, "DE_Results_GoodAnnotation.csv")
  
  rm(DE_Results_NoNA, DE_Results)
  
  print("Outputted object: DE_Results_GoodAnnotation")
}


#Applying the function to the DE results object:
FilteringDEResults_GoodAnnotation(DE_Results)



#This code includes a function within a function.

#This is the inner function:
GetContrastIDsforResultSet<-function(NamesOfFoldChangeColumns){
  #I split apart the column names:
  ColumnNames_BrokenUp<-strsplit(NamesOfFoldChangeColumns, "_")
  #Put them in a matrix format
  MatrixOfColumnNames_BrokenUp<-matrix(unlist(ColumnNames_BrokenUp), ncol=3,byrow=T)
  #And then grab the contrast ids:
  ContrastIDs_inCurrentDF<-MatrixOfColumnNames_BrokenUp[,2]
  rm(ColumnNames_BrokenUp, MatrixOfColumnNames_BrokenUp)
  return(ContrastIDs_inCurrentDF)
}

#This is the outer function:

ExtractingDEResultsForContrasts<-function(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts){
  
  print("These are all of the columns in the differential expression results for our current result set:")
  
  print(colnames(DE_Results_GoodAnnotation))
  # [1] "Probe"                       "NCBIid"                     
  # [3] "GeneSymbol"                  "GeneName"                   
  # [5] "pvalue"                      "corrected_pvalue"           
  # [7] "rank"                        "contrast_151617_coefficient"
  # [9] "contrast_151617_log2fc"      "contrast_151617_tstat"      
  # [11] "contrast_151617_pvalue"      "contrast_151618_coefficient"
  # [13] "contrast_151618_log2fc"      "contrast_151618_tstat"      
  # [15] "contrast_151618_pvalue"      "contrast_151619_coefficient"
  # [17] "contrast_151619_log2fc"      "contrast_151619_tstat"      
  # [19] "contrast_151619_pvalue" 
  
  print("These are the names of the Log(2) Fold Change Columns for our statistical contrasts of interest within the differential expression results for this particular result set:")
  
  NamesOfFoldChangeColumns<<-colnames(DE_Results_GoodAnnotation)[colnames(DE_Results_GoodAnnotation)%in%Contrasts_Log2FC]
  
  print(NamesOfFoldChangeColumns)
  #[1] "contrast_151617_log2fc" "contrast_151618_log2fc" "contrast_151619_log2fc"
  
  print("These are the names of the T-statistic Columns for our statistical contrasts of interest within the differential expression results for this particular result set:")
  
  NamesOfTstatColumns<<-colnames(DE_Results_GoodAnnotation)[colnames(DE_Results_GoodAnnotation)%in%Contrasts_Tstat]
  
  print(NamesOfTstatColumns)
  #[1] "contrast_151617_tstat" "contrast_151618_tstat" "contrast_151619_tstat"
  
  #Next we're going to pull out the contrast IDs associated with each result:
  
  ContrastIDs_inCurrentDF<-GetContrastIDsforResultSet(NamesOfFoldChangeColumns)
  
  print("These are the contrast ids for the statistical contrasts of interest within your current result set:")
  print(ContrastIDs_inCurrentDF)
  #[1] "151617" "151618" "151619"
  
  #Next we're going to grab some metadata to go with those statistical contrasts
  
  print("This is the dataset id for the result set and statistical contrasts:")
  Datasets_inCurrentDF<-ResultSet_contrasts$ExperimentID[ResultSet_contrasts$ContrastIDs%in%ContrastIDs_inCurrentDF]
  
  GSE_ID<<-Datasets_inCurrentDF[1]
  
  print(GSE_ID)
  #[1] "GSE126678"
  
  #And I would like the experimental factor information for our statistical contrast:
  Factors_inCurrentDF<-ResultSet_contrasts$ExperimentalFactors[ResultSet_contrasts$ContrastIDs%in%ContrastIDs_inCurrentDF]
  
  #We can combine those to make an interpretable unique identifier for each statistical comparison:
  ComparisonsOfInterest<<-paste(Datasets_inCurrentDF, Factors_inCurrentDF, sep="_" )
  
  print("These are the current names for your statistical contrasts of interest - if they are unwieldy, you may want to change them")
  print(ComparisonsOfInterest)
  
  #cleaning up the workspace:
  rm(Datasets_inCurrentDF, Factors_inCurrentDF, ContrastIDs_inCurrentDF)
  
}



#Read in the function:

#You can input this function by running the code discussed above to create the function in your R environment. 

#Alternatively, you can download the script for the function from our Github site and save the file in your working directory:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_ExtractingDEResultsForContrasts.R



#And then source it from your working directory:
#source("Function_ExtractingDEResultsForContrasts.R")




#Example function usage:
ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

"GSE248186_chronic inflammation with stress and derived from patient with major depressive disorder; chronic inflammation with stress and derived from patient with major depressive disorder"

#Those names are super unwieldy. I'm going to rename them:
#Note: Make sure the name still includes the GEO dataset id number (GSE...)
ComparisonsOfInterest<-c("GSE126678_CSDS_vs_Control")



CollapsingDEResults_OneResultPerGene<-function(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns){
  
  print("Double check that the vectors containing the fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:")
  
  print("# of rows with unique NCBI IDs:")
  print(length(unique(DE_Results_GoodAnnotation$NCBIid)))
  
  print("# of rows with unique Gene Symbols:")
  print(length(unique(DE_Results_GoodAnnotation$GeneSymbol)))
  
  #This makes a folder named after the dataset ID to store your results
  dir.create(paste("./", GSE_ID, sep=""))
  
  #And then sets the working directory to be that folder:
  setwd(paste("./", GSE_ID, sep=""))
  
  #Calculating the average Log2FC and T-stat for each gene:
  
  #Creating an empty list to store our results from each of our statistical contrasts:
  DE_Results_GoodAnnotation_FoldChange_Average<-list()
  DE_Results_GoodAnnotation_Tstat_Average<-list()
  DE_Results_GoodAnnotation_SE_Average<-list()
  
  #For each of the columns containing fold change and t-stat information for our statistical contrasts of interest:
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    
    #Grab our Log2FC column of interest:
    FoldChangeColumn<-select(DE_Results_GoodAnnotation, NamesOfFoldChangeColumns[i])
    
    #Grab our Tstat column of interest:
    TstatColumn<-select(DE_Results_GoodAnnotation, NamesOfTstatColumns[i])
    
    #Calculate the SE:
    DE_Results_GoodAnnotation_SE<-FoldChangeColumn[[1]]/TstatColumn[[1]]
    
    #Calculate the average Log2FC per gene:
    DE_Results_GoodAnnotation_FoldChange_Average[[i]]<-tapply(FoldChangeColumn[[1]], DE_Results_GoodAnnotation$NCBIid, mean)
    
    #Calculate the average Tstat per gene:
    DE_Results_GoodAnnotation_Tstat_Average[[i]]<-tapply(TstatColumn[[1]], DE_Results_GoodAnnotation$NCBIid, mean)
    
    #Calculate the average SE per gene:
    DE_Results_GoodAnnotation_SE_Average[[i]]<-tapply(DE_Results_GoodAnnotation_SE, DE_Results_GoodAnnotation$NCBIid, mean)
    
  }
  
  #Currently we have all of this averaged Log2FC information stored in a list, with one entry for each statistical contrast
  #Let's convert this to a data frame 
  DE_Results_GoodAnnotation_FoldChange_AveragedByGene<-do.call(cbind, DE_Results_GoodAnnotation_FoldChange_Average)
  
  print("Dimensions of Fold Change matrix, averaged by gene symbol:")
  print(dim(DE_Results_GoodAnnotation_FoldChange_AveragedByGene))
  
  #Name the columns in the dataframe in a manner that describes the dataset and factors for each statistic contrast
  colnames(DE_Results_GoodAnnotation_FoldChange_AveragedByGene)<-ComparisonsOfInterest
  
  #and then write it out to save it:
  write.csv(DE_Results_GoodAnnotation_FoldChange_AveragedByGene, "DE_Results_GoodAnnotation_FoldChange_AveragedByGene.csv")
  
  #And do the same for the T-stats:
  DE_Results_GoodAnnotation_Tstat_AveragedByGene<-do.call(cbind, DE_Results_GoodAnnotation_Tstat_Average)
  
  colnames(DE_Results_GoodAnnotation_Tstat_AveragedByGene)<-ComparisonsOfInterest
  
  write.csv(DE_Results_GoodAnnotation_Tstat_AveragedByGene, "DE_Results_GoodAnnotation_Tstat_AveragedByGene.csv")
  
  #And the same for the SE:
  DE_Results_GoodAnnotation_SE_AveragedByGene<-do.call(cbind, DE_Results_GoodAnnotation_SE_Average)
  
  colnames(DE_Results_GoodAnnotation_SE_AveragedByGene)<-ComparisonsOfInterest
  
  write.csv(DE_Results_GoodAnnotation_SE_AveragedByGene, "DE_Results_GoodAnnotation_SE_AveragedByGene.csv")
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  DE_Results_GoodAnnotation_SV<-(DE_Results_GoodAnnotation_SE_AveragedByGene)^2
  
  #Writing that information out to store it:
  write.csv(DE_Results_GoodAnnotation_SV, "DE_Results_GoodAnnotation_SV.csv")
  
  #Compiling all of our results into a single, big object 
  TempMasterResults<-list(Log2FC=DE_Results_GoodAnnotation_FoldChange_AveragedByGene, Tstat=DE_Results_GoodAnnotation_Tstat_AveragedByGene, SE=DE_Results_GoodAnnotation_SE_AveragedByGene, SV=DE_Results_GoodAnnotation_SV)
  
  #And then sending that object out into our environment outside the function
  #With a name that we can read and understand (DEResults for a particular GSEID)
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  #Letting us know what that name is:
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  #And then cleaning up our environment of temporary results:
  rm(TempMasterResults, DE_Results_GoodAnnotation, DE_Results_GoodAnnotation_SV, DE_Results_GoodAnnotation_SE, DE_Results_GoodAnnotation_FoldChange_AveragedByGene, DE_Results_GoodAnnotation_FoldChange_Average, DE_Results_GoodAnnotation_Tstat_AveragedByGene, DE_Results_GoodAnnotation_Tstat_Average, DE_Results_GoodAnnotation_SE_Average, FoldChangeColumn, TstatColumn, GSE_ID, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)
  
  #This sets the working directory back to the parent directory
  setwd("../")
  
}


#Reading in the function:

#You can input this function by running the code discussed above to create the function in your R environment. 

#Alternatively, you can download the script for the function from our Github site and save the file in your working directory:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_CollapsingDEResults_OneResultPerGene.R

library(dplyr)

#And then source it from your working directory:
source("Function_CollapsingDEResults_OneResultPerGene.R")

#Example usage:
CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)


##############################

DE_Results<-differentials[[2]]

#Applying the function to the DE results object:
FilteringDEResults_GoodAnnotation(DE_Results)

ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)

#Those names are super unwieldy. I'm going to rename them:
#Note: Make sure the name still includes the GEO dataset id number (GSE...)
ComparisonsOfInterest<-c("GSE150812_CSDSresilient_vs_Ctrl",  "GSE150812_CSDSsusceptible_vs_Ctrl")

CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)

##################

#Make sure you save your workspace and code!!! They will be separate files (one with .Rdata and one .R)
