
#A function for aligning all of our rat differential expression results from different datasets into a single data frame for Log2FCs and sampling variances (SVs):

AligningNewRatDatasets<-function(ListOfNewRatDEResults){
  
  #Making an empty list to hold our results:
  NewRat_MetaAnalysis_FoldChange_Dfs<-list()
  
  #Looping over all of the rat differential expression results:
  for(i in c(1:length(ListOfNewRatDEResults))){
    
    #Placing each of the log2FC results for each dataset into a single list
    #Each element in the list is formatted so that the rownames are Rat Entrez Gene ID and then there are columns containing the Log2FC for the differential expression results:
    NewRat_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Rat_GeneSymbol=row.names(ListOfNewRatDEResults[[i]][[1]]),ListOfNewRatDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  #Letting the user know the structure of the set of Log2FC dataframes that we are starting out with:
  print("NewRat_MetaAnalysis_FoldChange_Dfs:")
  print(str(NewRat_MetaAnalysis_FoldChange_Dfs))
  
  #Running "join all" on the list to align all of the results by Entrez Gene ID and make them a single data frame:
  NewRat_MetaAnalysis_FoldChanges<<-join_all(NewRat_MetaAnalysis_FoldChange_Dfs, by="Rat_GeneSymbol", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  #Letting the user know the structure of the dataframe that we just created:
  print("NewRat_MetaAnalysis_FoldChanges:")
  print(str(NewRat_MetaAnalysis_FoldChanges))
  
  #Doing the same steps for the sampling variances:
  
  NewRat_MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfNewRatDEResults))){
    NewRat_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Rat_GeneSymbol=row.names(ListOfNewRatDEResults[[i]][[4]]),ListOfNewRatDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("NewRat_MetaAnalysis_SV_Dfs:")
  print(str(NewRat_MetaAnalysis_SV_Dfs))
  
  NewRat_MetaAnalysis_SV<<-join_all(NewRat_MetaAnalysis_SV_Dfs, by="Rat_GeneSymbol", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("NewRat_MetaAnalysis_SV:")
  print(str(NewRat_MetaAnalysis_SV))
  
  #Cleaning up our environment to remove unneeded objects:
  rm(NewRat_MetaAnalysis_SV_Dfs, NewRat_MetaAnalysis_FoldChange_Dfs)
}


#Then join 
NewRat_MetaAnalysis_FoldChanges
#to Rat_MetaAnalysis_FoldChanges
#using Hom_MouseVsRat to translate between EntrezId and GeneSymbol
#We will probably need to rename column names in Hom_MouseVsRat to do this to make them match our data.frame

#And do the same for this...
NewRat_MetaAnalysis_SV
