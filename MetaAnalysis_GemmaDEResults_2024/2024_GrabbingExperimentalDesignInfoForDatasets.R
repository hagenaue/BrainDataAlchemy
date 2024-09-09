#Extracting Experimental Design details from Gemma's database
#Megan Hagenauer, July 8 2024

#The Gemma database includes differential expression results for each of the independent variables included in the transcriptional profiling datasets

#We will need to review the design of these differential expression analyses for two reasons:

##############################

#First Goal:

#Sometimes it is difficult to determine from the abstract and dataset annotation which variables were *actually manipulated* as part of the experiment. 

#This can be especially true for datasets that were part of much larger studies containing multiple experiments

#Gemma's information about the experimental design for the transcriptional profiling experiment can clarify this. 

#This code pulls out that information and places it in a flat format (dataset=row, column=contrast metadata) that can be easily combined with the spreadsheet file where you are already keeping your dataset inclusion/exclusion notes 

#First, you will need to read in a vector of your dataset ids 
#You can either write them out by hand, e.g., 
#ExperimentIDs<-c("GSE135306", "GSE81672", "GSE180465", "GSE179667", "GSE146358", "GSE128255", "GSE187418")

#... or you can copy them into a spreadsheet as a single column and save them as a .csv file (preferably in the same order as your inclusion/exclusion spreadsheet, without a column name (header)) and read them in:
ExperimentIDs<-read.csv("ExampleExperimentIDs.csv", header=FALSE, stringsAsFactors = FALSE)
str(ExperimentIDs)
#This is a dataframe.
#... but the following code needs a vector, so let's grab the first (and only) column
#And make it a vector
ExperimentIDs<-ExperimentIDs[,1]
str(ExperimentIDs)

#Next, we'll make some empty vectors to store the information that we are going to collect about each dataset:
FactorInfo<-vector(mode="character", length=length(ExperimentIDs))
BaselineFactorInfo<-vector(mode="character", length=length(ExperimentIDs))
TreatmentFactorInfo<-vector(mode="character", length=length(ExperimentIDs))

#And combine them into a nice dataframe:
FactorInfoDF<-cbind.data.frame(ExperimentIDs, FactorInfo, TreatmentFactorInfo, BaselineFactorInfo)

#We'll then loop over each of the datasets:

for(i in c(1:length(ExperimentIDs))){
  
  #For each dataset, we will use Gemma's API to access the experimental design info:
  Design<-gemma.R::get_dataset_differential_expression_analyses(ExperimentIDs[i])
  
  if(nrow(Design)>0){
    
  #For the dataset, we'll grab the vector of factor categories included in the design
  #And collapse it to a single entry in our data frame.
  FactorInfoDF[i,2]<-paste(Design$factor.category, collapse="; ")
  
  #For the dataset, we'll next grab the factors included in each result set for the dataset 
  #To do this, we'll first make a vector to store information
  ExperimentalFactorVector<-vector(mode="character", length=1)
  
  #And then loop over each of the result sets for the dataset:
  for(j in c(1:length(Design$result.ID))){
    
    #And grab the factors for the result set and concatenate them to our vector of all factors included in any result set for the dataset:
    ExperimentalFactorVector<-c(ExperimentalFactorVector,Design$experimental.factors[[j]]$summary)
  }
  
  #Then we'll remove the first (empty) entry in the vector and collapse the vector of factors into a single entry for our data frame:
  FactorInfoDF[i,3]<-paste(ExperimentalFactorVector[-1], collapse="; ")
  
  #For the dataset, we'll next grab the baseline/control/reference for each of the factors included in each result set for the dataset 
  #To do this, we'll first make a vector to store information
  BaselineFactorVector<-vector(mode="character", length=1)
  
  #And then loop over each of the result sets for the dataset:
  for(k in c(1:length(Design$result.ID))){
    
    #And grab the baseline/control/reference definitions for the result set and concatenate them to our vector of all baselines included in any result set for the dataset:
    BaselineFactorVector<-c(BaselineFactorVector,Design$baseline.factors[[k]]$summary)
  }
  
  #Then we'll remove the first (empty) entry in the vector and collapse the vector of baselines into a single entry for our data frame:
  FactorInfoDF[i,4]<-paste(BaselineFactorVector[-1], collapse="; ")
  
  #And clean up our workspace before we start the loop again
  rm(Design, ExperimentalFactorVector, BaselineFactorVector)
  
  }else{
    rm(Design)
  }
}

#You can write out this object and use it to help screen datasets
write.csv(FactorInfoDF, "FactorInfoDF.csv")

#################################

#Second Goal:

#For the meta-analysis, we will be extracting the differential expression results from Gemma. The differential expression results for each dataset may include multiple statistical contrasts (e.g., drug1 vs. vehicle, drug2 vs. vehicle). 
#Each of these contrasts are labeled with a result id and contrast id within the Gemma database. 
#We will need to know which of these ids are relevant to our project goals to easily extract their results.

#We will also need to double-check that these statistical contrasts are set up in a manner that makes sense for our experiments:

#First, for experiments that include more than one brain region, we will need to double-check that the results have been subsetted by brain region (instead of including brain region ("OrganismPart") as a factor in the model). If they haven't been subsetted by region, we will probably need to re-run the differential expression analysis.

#Depending on the goals of the meta-analysis, we may also need to re-run the differential expression analysis to remove other unwanted subjects (e.g., removing subjects with genotypes that might interfere with our results)

#Second, we will need to double-check that the comparisons include an appropriate reference group - sometimes they are reversed in Gemma (e.g., having the drug treatment set as the baseline, with vehicle as the manipulation). If this is the case, we will need to invert the effects when we input them into our meta-analysis (multiply the effects by -1).


#This code makes a data.frame that includes all of the contrast ids for each dataset with their basic metadata in a format that is easily readable in a spreadsheet program.

#Making an empty data.frame to store results:

ResultSets_toScreen<-data.frame(ExperimentID="NA",ResultSetIDs="NA", ContrastIDs="NA", ExperimentIDs="NA", FactorCategory="NA", ExperimentalFactors="NA", BaselineFactors="NA", Subsetted=FALSE, SubsetBy="NA")

str(ResultSets_toScreen)
# 'data.frame':	1 obs. of  9 variables:
# $ ExperimentID       : chr "NA"
# $ ResultSetIDs       : chr "NA"
# $ ContrastIDs        : chr "NA"
# $ ExperimentIDs      : chr "NA"
# $ FactorCategory     : chr "NA"
# $ ExperimentalFactors: chr "NA"
# $ BaselineFactors    : chr "NA"
# $ Subsetted          : logi FALSE
# $ SubsetBy           : chr "NA"

#We will then loop over each of the datasets:

for(i in c(1:length(ExperimentIDs))){
  
  #For each dataset, we will use Gemma's API to access the experimental design info:
  Design<-gemma.R::get_dataset_differential_expression_analyses(ExperimentIDs[i])

  if(nrow(Design)>0){
#Next, we'll make some empty vectors to store the experimental factor and baseline factor information for each result id for the dataset:
ExperimentalFactors<-vector(mode="character", length(Design$result.ID))
BaselineFactors<-vector(mode="character", length(Design$result.ID))

#We will then loop over each of the result ids for the dataset:
for(j in c(1:length(Design$result.ID))){
  
  #And grab the vector of experimental factors associated with that result id
  ExperimentalFactorVector<-Design$experimental.factors[[j]]$summary
  #And collapse that info down to a single entry that will fit in our data.frame
  ExperimentalFactors[j]<-paste(ExperimentalFactorVector, collapse="; ")

  #And then grab the vector of baseline/control/reference values associated with that result id
  BaselineFactorVector<-Design$baseline.factors[[j]]$summary
  #And collapse that info down to a single entry that will fit in our data.frame
  BaselineFactors[j]<-paste(BaselineFactorVector, collapse="; ")
}

#Some of the datasets are subsetted for the differential expression analyses
#We will make an empty vector to store subset information for each result id
SubsetBy<-vector(mode="character", length(Design$result.ID))

#Then we will determine whether the dataset is subsetted:
if(Design$isSubset[1]==TRUE){
  
  #If it is subsetted, we will loop over each result id for the dataset
  for (j in c(1:length(Design$result.ID))){
    
    #And grab the vector of subsetting information
    SubsetByVector<-Design$subsetFactor[[j]]$summary
    
    #And then collapse that information down to a single entry that will fit in our dataframe
    SubsetBy[j]<-paste(SubsetByVector, collapse="; ")
  }  
  
  #if the dataset wasn't subsetted for the differential expression analysis:
}else{
  #We'll just make a vector of NA values to put in the "Subsetted by" column
  SubsetBy<-rep(NA, length((Design$result.ID)))
}

#Then we combine all of the information for all of the result sets for the dataset into a dataframe
ResultSets_ForExperiment<-cbind.data.frame(ExperimentID=rep(ExperimentIDs[i],length(Design$result.ID)),ResultSetIDs=Design$result.ID, ContrastIDs=Design$contrast.ID, ExperimentIDs=Design$experiment.ID, FactorCategory=Design$factor.category, ExperimentalFactors, BaselineFactors, Subsetted=Design$isSubset, SubsetBy)

#And add that information as rows to our data frame including the result set information for all datasets:
ResultSets_toScreen<-rbind.data.frame(ResultSets_toScreen, ResultSets_ForExperiment)

#Then clean up our space before looping to the next dataset:
rm(ResultSets_ForExperiment, Design, ExperimentalFactors, BaselineFactors, SubsetBy)
  
  }else{
  rm(Design)
}
  
}

#When we're done, we'll want to remove the initial (empty) row in our data.frame:
ResultSets_toScreen<-ResultSets_toScreen[-1,]

#We can make some empty vectors that we can use to store screening notes:
Include<-vector(mode="character", length=nrow(ResultSets_toScreen))
WrongBaseline<-vector(mode="character", length=nrow(ResultSets_toScreen))
ResultsNotRegionSpecific<-vector(mode="character", length=nrow(ResultSets_toScreen))
ReAnalyze<-vector(mode="character", length=nrow(ResultSets_toScreen))

#And add them as columns to our dataframe:            
ResultSets_toScreen<-cbind.data.frame(ResultSets_toScreen, Include, WrongBaseline, ResultsNotRegionSpecific, ReAnalyze)

#And then write everything out as a .csv file that we can easily mark up in a spreadsheet program:
write.csv(ResultSets_toScreen, "ResultSets_toScreen.csv")


