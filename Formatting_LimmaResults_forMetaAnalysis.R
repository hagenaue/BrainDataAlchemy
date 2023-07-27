#Preparing our limma results for meta-analysis
#Megan Hagenauer
#July 27 2023

###############

#Start with an empty workspace

###############

#Cleaning up the annotation for our results

#Example from Anne's results
setwd("/Users/hagenaue/Documents/drive-download-20230727T135146Z-001")

TempResultsJoined<-read.csv("Annotation_LimmaResults_GSE92718_all.csv", header=TRUE, stringsAsFactors = FALSE)
str(TempResultsJoined)

#################

#Reading in the function

FilteringDEResults_GoodAnnotation<-function(TempResultsJoined){
  
  print("# of rows in results")
  print(nrow(TempResultsJoined))
  
  print("# of rows with missing NCBI annotation:")
  print(sum(TempResultsJoined$NCBIid==""|TempResultsJoined$NCBIid=="null"))
  
  print("# of rows with NA NCBI annotation:")
  print(sum(is.na(TempResultsJoined$NCBIid)))
  
  print("# of rows with missing Gene Symbol annotation:")
  print(sum(TempResultsJoined$GeneSymbol==""|TempResultsJoined$GeneSymbol=="null"))
  
  print("# of rows mapped to multiple NCBI_IDs:")
  print(length(grep('\\|', TempResultsJoined$NCBIid)))
  
  print("# of rows mapped to multiple Gene Symbols:")
  print(length(grep('\\|', TempResultsJoined$GeneSymbol)))
  
  #I only want the subset of data which contains rows that do not contain an NCBI EntrezID of ""
  TempResultsJoined_NoNA<-TempResultsJoined[(TempResultsJoined$NCBIid==""|TempResultsJoined$NCBIid=="null")==FALSE & is.na(TempResultsJoined$NCBIid)==FALSE,]
  
  if(length(grep('\\|', TempResultsJoined_NoNA$NCBIid))==0){
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA
  }else{
    #I only want rows annotated with a single Gene Symbol (no pipe):
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA[-(grep('\\|', TempResultsJoined_NoNA$NCBIid)),]
  }
  
  print("# of rows with good annotation")
  print(nrow(TempResultsJoined_NoNA_NoMultimapped))
  
  #For record keeping (sometimes useful for troubleshooting later)
  write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_NoNA_NoMultimapped.csv")
  
  rm(TempResultsJoined_NoNA)
  
  print("Outputted object: TempResultsJoined_NoNA_NoMultimapped")
}

#################

#Example of using the function for a dataset

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] NA
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 5
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 16867
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"


#################

#Adapting the Log2FC and SE extraction from last year:

#We need to make a matrix of the Log2FoldChange values for the comparisons of interest. 
#These columns start with "Coef" and are specific to your comparisons of interest


colnames(TempResultsJoined_NoNA_NoMultimapped)

TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.SumExp_Subset_noBad_GSE92718_Filtered.treatment_factorcuff.operation)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#We also need to make a matrix of the Tstat values for the comparisons of interest
#These columns start with "t." and are specific to your comparisons of interest

TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.SumExp_Subset_noBad_GSE92718_Filtered.treatment_factorcuff.operation)

str(TempResultsJoined_NoNA_NoMultimapped_Tstats)

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid


#Let's rename our columns to something nicer

ComparisonsOfInterest<-c("CuffOperation_vs_Ctrl")

colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest

#Sanity check:
#1.386/8.37



###################

#Reading in the Extracting Results function:

ExtractingDEResults<-function(GSE_ID, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats){
  
  #Let's rename our columns to something nicer
  
  #We calculate the standard error by dividing the log2FC by the tstat
  TempResultsJoined_NoNA_NoMultimapped_SE<-TempResultsJoined_NoNA_NoMultimapped_FoldChanges/TempResultsJoined_NoNA_NoMultimapped_Tstats
  str(TempResultsJoined_NoNA_NoMultimapped_SE)
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  TempResultsJoined_NoNA_NoMultimapped_SV<-(TempResultsJoined_NoNA_NoMultimapped_SE)^2
  str(TempResultsJoined_NoNA_NoMultimapped_SV)
  
  TempMasterResults<-list(Log2FC=TempResultsJoined_NoNA_NoMultimapped_FoldChanges, Tstat=TempResultsJoined_NoNA_NoMultimapped_Tstats, SE=TempResultsJoined_NoNA_NoMultimapped_SE, SV=TempResultsJoined_NoNA_NoMultimapped_SV)
  
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  rm(TempMasterResults, TempResultsJoined_NoNA_NoMultimapped_SV, TempResultsJoined_NoNA_NoMultimapped_SE, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)
  
}


###################


#Example usage for the extracting results function:

ExtractingDEResults(GSE_ID="GSE92718", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)
# num [1:16862, 1] 0.1655 0.0943 0.1069 0.1443 0.0382 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16862] "13654" "215418" "13656" "19252" ...
# ..$ : chr "CuffOperation_vs_Ctrl"
# num [1:16862, 1] 0.02739 0.00889 0.01142 0.02083 0.00146 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:16862] "13654" "215418" "13656" "19252" ...
# ..$ : chr "CuffOperation_vs_Ctrl"
# [1] "Output: Named DEResults_GSE92718"

str(DEResults_GSE92718)
# List of 4
# $ Log2FC: num [1:16862, 1] 1.386 0.668 0.72 0.974 0.245 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:16862] "13654" "215418" "13656" "19252" ...
# .. ..$ : chr "CuffOperation_vs_Ctrl"
# $ Tstat : num [1:16862, 1] 8.37 7.09 6.74 6.75 6.41 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:16862] "13654" "215418" "13656" "19252" ...
# .. ..$ : chr "CuffOperation_vs_Ctrl"
# $ SE    : num [1:16862, 1] 0.1655 0.0943 0.1069 0.1443 0.0382 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:16862] "13654" "215418" "13656" "19252" ...
# .. ..$ : chr "CuffOperation_vs_Ctrl"
# $ SV    : num [1:16862, 1] 0.02739 0.00889 0.01142 0.02083 0.00146 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:16862] "13654" "215418" "13656" "19252" ...
# .. ..$ : chr "CuffOperation_vs_Ctrl"

#Sanity Check:

#Calculating the SE from the Log2FC and Tstat for the first gene(row)
1.386/8.37
#[1] 0.1655914

#Calculating the SV for the first gene
0.1655914^2
#[1] 0.02742051
#Looks good (the decimal differences are probably just due to rounding)

#Calculating the SE from the Log2FC and Tstat for the second gene(row)
0.668/7.09
#[1] 0.09421721

#Calculating the SV for the second gene
0.09421721^2
#[1] 0.008876883
#Looks good (the decimal differences are probably just due to rounding)


#################

#Cleaning up our environment before moving on to the next dataset:

rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)

#################

#Running the next dataset:

##

#This part of the code will be dataset specific:

#Setting the working directory to where the limma results for the dataset are located:
setwd("/Users/hagenaue/Documents/drive-download-20230727T135232Z-001")

#Depending on file format, you may need to use read.delim instead of read.csv:
TempResultsJoined<-read.delim("Limma_interactionresults_85136.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(TempResultsJoined)

###

#This part is not dataset specific:
FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows in results"
# [1] 20671
# [1] "# of rows with missing NCBI annotation:"
# [1] 14
# [1] "# of rows with NA NCBI annotation:"
# [1] 0
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 14
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 7
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 7
# [1] "# of rows with good annotation"
# [1] 20650
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

colnames(TempResultsJoined_NoNA_NoMultimapped)

##

#This part is (somewhat) dataset specific
#We just need to replace the column names for the Log2FC ("Coef") and T-stats ("t.")
#With the appropriate names for this dataset:

TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.interaction.stress.x.female)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#We also need to make a matrix of the Tstat values for the comparisons of interest
#These columns start with "t." and are specific to your comparisons of interest

TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.interaction.stress.x.female)

str(TempResultsJoined_NoNA_NoMultimapped_Tstats)

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:

ComparisonsOfInterest<-c("StressEffects_InteractionWSex")

colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest

#This is dataset specific, only because we need to provide the GSE#:

ExtractingDEResults(GSE_ID="GSE85136", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)

##

#This part is not dataset specific:

str(DEResults_GSE85136)

rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)

##################


#Move on to the next dataset...