
GettingResultSetInfoForDatasets<-function(ExperimentIDs){
  
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
  
  print("The Result Sets for your Datasets have been outputted into ResultSets_toScreen.csv")
  print(str(ResultSets_toScreen))
  
  #And clean up our environment:
  rm(Include, WrongBaseline, ResultsNotRegionSpecific, ReAnalyze)
}