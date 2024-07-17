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