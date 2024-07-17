
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
