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