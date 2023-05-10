#Annaka Saffron, summer 2022

#Code for identifying useful datasets:

# Search Terms
# LPS
# 
# Bacterial infection 
# 
# LPS hippocampus 
# R code

if(!requireNamespace("devtools", quietly=T)){
  install.packages("devtools")
}

devtools:: install_github("PavlidisLab/gemma.R", force=T)


setwd("C:/Users/annak/Documents/gemma datasets")


lps1 <- gemma.R :: searchDatasets("LPS", taxon= "rat", limit= 100)
write.table(x=lps1, file= "lps search result1.txt", sep= "\t")

##11 rows 

lps2 <- gemma.R :: searchDatasets("LPS", taxon= "mouse", limit= 100)
write.table(x=lps2, file= "lps search result2.txt", sep= "\t")

lps3 <- gemma.R :: searchDatasets("LPS", taxon= "mouse", offset=100,  limit= 100)
write.table(x=lps3, file= "lps search result3.txt", sep= "\t")


lps4 <- gemma.R :: searchDatasets("LPS", taxon= "mouse", offset=200,  limit= 100)
write.table(x=lps4, file= "lps search result4.txt", sep= "\t")
##229 rows for lps mouse
##had to offset (skip specified num of objects found from database) to avoid the
##100 obj limit ex.) lps2 is the first 100 mouse lps studies, lps3 is the 101-200
##lps mouse studies and lps4 is the 201-229 mouse lps studies
##none of the other search terms hit the limit so they only have 1 and 2 (rat, mouse)

lps_bacteria1 <- gemma.R :: searchDatasets("bacterial infection", taxon= "rat", limit= 100)
write.table(x=lps_bacteria1, file= "bacteria search result1.txt", sep= "\t")
##1 row
lps_bacteria2 <- gemma.R :: searchDatasets("bacterial infection", taxon="mouse", limit= 100)
write.table(x=lps_bacteria2, file= "bacteria search result2.txt", sep= "\t")
##71 row

lps_hipp1 <- gemma.R :: searchDatasets("lps hippocampus", taxon="rat", limit=100)
write.table(x=lps_hipp1, file= "lpshippocampussearch result1.txt", sep= "\t")
##0 rows
lps_hipp2 <- gemma.R :: searchDatasets("lps hippocampus", taxon="mouse", limit=100)
write.table(x=lps_hipp2, file= "lpshippocampussearch result2.txt", sep= "\t")
##14 rows

combined<-rbind.data.frame(lps1, lps2, lps3, lps4, lps_hipp1, lps_hipp2, lps_bacteria1, lps_bacteria2)

##326 rows
unique_combined <- unique(combined) ##304 rows

write.table(x=unique_combined, file= "unique_data.txt", sep= "\t")


##date= 7/11/2022

#*then used excel conditional formatting and ctrl F to pull only RNA seq or microarray data from hippocampus (including ammon’s horn, DG, CA1,2,3,4)* this resulted in only 5 datasets.

#################

#Code for reading in datasets and meta-analysis:

##packages


library(plyr)
library(metafor)

library(reshape)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest", force=TRUE)
library(multtest)

library(MAd)
##Functions for reading in




ReadingInGemmaDE<-function(ResultSetFileNames){
  
  #Reading in the analysis results file:
  TempAnalysisResults<-read.delim("analysis.results.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")
  
  #I placed these results in a list format:
  TempResultsToJoin<-list(TempAnalysisResults)
  
  for(i in c(1:length(ResultSetFileNames))){
    #The result set files are the detailed differential expression results for a particular variable:
    TempResultsToJoin[[i]]<-read.delim(ResultSetFileNames[i], sep="\t", stringsAsFactors = FALSE, comment.char = "#")
  }
  
  TempResultsJoined<<-join_all(TempResultsToJoin, by="Element_Name")
  #Note: I've heard from other students that the join() and join_all() functions in the plyr package can cause problems in newer versions of R - you may need to replace this with merge and merge_all
  
  #Saving the joined results:
  write.csv(TempResultsJoined, "TempResultsJoined.csv")
  
  rm(TempAnalysisResults, TempResultsToJoin)
  
  print("Outputted object: TempResultsJoined")
}


FilteringDEResults_GoodAnnotation<-function(TempResultsJoined){
  
  print("# of rows with missing NCBI annotation:")
  print(sum(TempResultsJoined$NCBI_ID==""|TempResultsJoined$NCBI_ID=="null"))
  
  print("# of rows with missing Gene Symbol annotation:")
  print(sum(TempResultsJoined$Gene_Symbol==""|TempResultsJoined$Gene_Symbol=="null"))
  
  print("# of rows mapped to multiple NCBI_IDs:")
  print(length(grep('\\|', TempResultsJoined$NCBI_ID)))
  
  print("# of rows mapped to multiple Gene Symbols:")
  print(length(grep('\\|', TempResultsJoined$Gene_Symbol)))
  
  #I only want the subset of data which contains rows that do not contain a Gene Symbol of ""
  TempResultsJoined_NoNA<-TempResultsJoined[(TempResultsJoined$Gene_Symbol==""|TempResultsJoined$Gene_Symbol=="null")==FALSE,]
  
  if(length(grep('\\|', TempResultsJoined_NoNA$Gene_Symbol))==0){
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA
  }else{
    #I only want rows annotated with a single Gene Symbol (no pipe):
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA[-(grep('\\|', TempResultsJoined_NoNA$Gene_Symbol)),]
  } 
  print("# of rows with good annotation")
  print(nrow(TempResultsJoined_NoNA_NoMultimapped))
  
  #For record keeping (sometimes useful for troubleshooting later)
  write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_NoNA_NoMultimapped.csv")
  
  rm(TempResultsJoined_NoNA, TempResultsJoined_NoNA_NoMultimapped)
  
  print("Outputted object: TempResultsJoined_NoNA_NoMultimapped")
}



CollapsingDEResults_OneResultPerGene<-function(GSE_ID, TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns){
  
  print("Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:")
  
  print("# of rows with unique NCBI IDs:")
  print(length(unique(TempResultsJoined_NoNA_NoMultimapped$NCBI_ID)))
  
  print("# of rows with unique Gene Symbols:")
  print(length(unique(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol)))
  
  #We will need both the Log2FC and T-stats averaged:
  
  TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average<-list()
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    
    TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average[[i]]<-tapply(NamesOfFoldChangeColumns[i][[1]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    
  }
  
  TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average)
  
  print("Dimensions of Fold Change matrix, averaged by gene symbol:")
  print(dim(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol))
  
  colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol)<-ComparisonsOfInterest
  
  write.csv(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol.csv")
  
  
  TempResultsJoined_NoNA_NoMultimapped_Tstat_Average<-list()
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    
    TempResultsJoined_NoNA_NoMultimapped_Tstat_Average[[i]]<-tapply(NamesOfTstatColumns[i][[1]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    
  }
  
  TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_Tstat_Average)
  
  colnames(TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol)<-ComparisonsOfInterest
  
  write.csv(TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol.csv")
  
  
  
  TempResultsJoined_NoNA_NoMultimapped_SE<-list()
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    TempResultsJoined_NoNA_NoMultimapped_SE[[i]]<-NamesOfFoldChangeColumns[i][[1]]/NamesOfTstatColumns[i][[1]]
  }
  
  
  TempResultsJoined_NoNA_NoMultimapped_SE_Average<-list()
  
  for(i in c(1:length(NamesOfFoldChangeColumns))){
    
    TempResultsJoined_NoNA_NoMultimapped_SE_Average[[i]]<-tapply(TempResultsJoined_NoNA_NoMultimapped_SE[[i]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    
  }
  
  TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_SE_Average)
  
  colnames(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol)<-ComparisonsOfInterest
  
  write.csv(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol.csv")
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  TempResultsJoined_NoNA_NoMultimapped_SV<-(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol)^2
  
  write.csv(TempResultsJoined_NoNA_NoMultimapped_SV, "TempResultsJoined_NoNA_NoMultimapped_SV.csv")
  
  TempMasterResults<-list(Log2FC=TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, Tstat=TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, SE=TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol, SV=TempResultsJoined_NoNA_NoMultimapped_SV)
  
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  rm(TempMasterResults, TempResultsJoined_NoNA_NoMultimapped_SV, TempResultsJoined_NoNA_NoMultimapped_SE, TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average, TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_Tstat_Average)
  
}




AligningDatasets<-function(ListOfDEResults){
  
  MetaAnalysis_FoldChange_Dfs<-list()
  
  for(i in c(1:length(ListOfDEResults))){
    MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[1]]),ListOfDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  print("MetaAnalysis_FoldChange_Dfs:")
  print(str(MetaAnalysis_FoldChange_Dfs))
  
  MetaAnalysis_FoldChanges<<-join_all(MetaAnalysis_FoldChange_Dfs, by="x", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("MetaAnalysis_FoldChanges:")
  print(str(MetaAnalysis_FoldChanges))
  
  MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfDEResults))){
    MetaAnalysis_SV_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[4]]),ListOfDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("MetaAnalysis_SV_Dfs:")
  print(str(MetaAnalysis_SV_Dfs))
  
  MetaAnalysis_SV<<-join_all(MetaAnalysis_SV_Dfs, by="x", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("MetaAnalysis_SV:")
  print(str(MetaAnalysis_SV))
  
  rm(MetaAnalysis_SV_Dfs, MetaAnalysis_FoldChange_Dfs)
}





RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 6)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
  for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
    
    if(skip_to_next){}else{
      TempMeta<-rma(effect, var)
      metaOutput[i, 1]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, 2]<-TempMeta$se #gives standard error
      metaOutput[i, 3]<-TempMeta$pval #gives pval
      metaOutput[i, 4]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, 5]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 6]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons")
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta$x
  
  metaOutput<<-metaOutput
  return(metaOutput)
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}


FalseDiscoveryCorrection<-function(metaOutput){
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[7]<-"FDR"
  
  metaOutputFDR<<-metaOutputFDR
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  write.csv(metaOutputFDR, "metaOutputFDR.csv")
  
  #a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR[order(metaOutputFDR[,3]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval_wHDRFData.csv")
  
  print("Do we have any genes that are statistically significant following false discovery rate correction?")
  print(sum(metaOutputFDR[,7]<0.10, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR[order(metaOutputFDR[,3]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
}





MakeForestPlots<-function(GeneSymbol){
  
  pdf(paste("ForestPlot_", GeneSymbol, ".pdf", sep=""), height=5, width=8)
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  
  forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-1],  xlim=c(-3, 3))
  
  mtext(paste(GeneSymbol), line=-1.5, cex=2)
  dev.off()
}




RunMetaAnalysisDiagnostics<-function(GeneSymbol){
  
  print(GeneSymbol)  
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  var<-as.numeric(MetAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  
  FAILED <- FALSE
  tryCatch(TempMeta<-rma(effect, var), error = function(e) {FAILED <<- TRUE})
  
  if(FAILED){print("Meta-Analysis Failed")}else{
    TempMeta<-rma(effect, var)
    print(TempMeta)
    print("# of Comparisons:")
    print(8-sum(is.na(effect)))
    print("Cook's Distance:")
    print(data.frame(Study=colnames(MetaAnalysis_FoldChanges_ForMeta )[-1][is.na(effect)==FALSE],Distance=cooks.distance(TempMeta)))
    print("DF Betas:")
    print(data.frame(Study=colnames(MetaAnalysis_FoldChanges_ForMeta )[-1][is.na(effect)==FALSE], DfBetas=dfbetas(TempMeta)))
    rm(TempMeta)
  }
  rm(effect, var)
}



DistanceMetrics<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 3)
CooksDistanceForStudies<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 8)
#And then run a loop that run's a meta-analysis on the differential expression results (columns 2-7) for each gene (row):
for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
  var<-as.numeric(MetAnalysis_SV_ForMeta[i,-1])
  
  skip_to_next <- FALSE
  tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
  
  if(skip_to_next){}else{
    TempMeta<-rma(effect, var)
    DistanceMetrics[i, 1]<-TempMeta$QEp  
    DistanceMetrics[i, 2]<-TempMeta$I2 
    CooksDistance<-cooks.distance(TempMeta)
    DistanceMetrics[i, 3]<-max(CooksDistance) 
    CooksDistanceForStudies[i, as.numeric(names(cooks.distance(TempMeta)))]<-CooksDistance
    rm(TempMeta, CooksDistance)
  }
  rm(effect, var)
}



##gse12284
getwd()
setwd("C:/Users/annak/Documents/R/GSE12284")
list.files()

ReadingInGemmaDE(ResultSetFileNames="resultset_ID471315.data.txt")

str(TempResultsJoined)

FilteringDEResults_GoodAnnotation(TempResultsJoined)

str(TempResultsJoined_NoNA_NoMultimapped)

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE12284", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest = c("LPS_vs_Ctrl12284"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_lipopolysaccharide), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_lipopolysaccharide))


##GSE126678

getwd()
setwd("C:/Users/annak/Documents/R/GSE126678/15127_GSE126678_diffExpAnalysis_119797")
list.files()

ReadingInGemmaDE(ResultSetFileNames="resultset_ID488449.data.txt" )

str(TempResultsJoined)

FilteringDEResults_GoodAnnotation(TempResultsJoined)

str(TempResultsJoined_NoNA_NoMultimapped)

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE1266781", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest = c("LPSAcute_vs_Ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_lipopolysaccharide_Acute.immune.challenge._vehicle), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_lipopolysaccharide_Acute.immune.challenge._vehicle))
##for acute v s vehicle
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE1266782", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest = c("LPSAcute_vs_LPSChr"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_Long.term.subchronic.immune.challenge.acute.immune.challenge_lipopolysaccharide_lipopolysaccharide), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_Long.term.subchronic.immune.challenge.acute.immune.challenge_lipopolysaccharide_lipopolysaccharide))
##for acute vs subchronic

##GSE181285
getwd()
setwd("C:/Users/annak/Documents/R/GSE181285/21341_GSE181285_diffExpAnalysis_216429")
list.files()

ReadingInGemmaDE(ResultSetFileNames="resultset_ID507832.data.txt" )

str(TempResultsJoined)

FilteringDEResults_GoodAnnotation(TempResultsJoined)

str(TempResultsJoined_NoNA_NoMultimapped)

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE181285", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest = c("LPS_vs_Ctrl181285"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_850.g.kg_lipopolysaccharide), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_850.g.kg_lipopolysaccharide))


##81024

getwd()

setwd("13018_GSE81024_diffExpAnalysis_99567")
list.files()
#This procedure for reading in the differential expression results was pretty easily generalizable, so I functionalized it:

ReadingInGemmaDE("resultset_ID480667.data.txt")

setwd("C:/Users/annak/Documents/R/GSE81024")

MoreAnnotation<-read.delim("GPL7202-9760.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
colnames(MoreAnnotation)
NeededAnnotation<-MoreAnnotation[,c(1,7:8)]
colnames(TempResultsJoined)
colnames(NeededAnnotation)<-c("Element_Name", "Gene_Symbol_New","Gene_Name_New")
library(plyr)
TempResultsJoined_wAnnotation<-join(TempResultsJoined, NeededAnnotation, by="Element_Name", type="left")
TempResultsJoined_wAnnotation$Gene_Symbol<-TempResultsJoined_wAnnotation$Gene_Symbol_New
TempResultsJoined_wAnnotation$Gene_Name<-TempResultsJoined_wAnnotation$Gene_Name_New
TempResultsJoined<-TempResultsJoined_wAnnotation

FilteringDEResults_GoodAnnotation(TempResultsJoined)



CollapsingDEResults_OneResultPerGene(GSE_ID="GSE81024", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest = c("LPS_vs_Ctrl81024"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_lipopolysaccharide), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_lipopolysaccharide))


rm(TempResultsJoined)
rm(TempResultsJoined_NoNA_NoMultimapped)
rm(TempResultsJoined_wAnnotation)









ListOfDEResults<-list(DEResults_GSE12284, DEResults_GSE1266781, DEResults_GSE1266782, DEResults_GSE181285, DEResults_GSE81024)
#I just discovered that if I have datasets with *comparisons with the same name*, when I join the datasets some of those comparisons disappear. :(

AligningDatasets(ListOfDEResults)


NumberOfComparisons=5
CutOffForNAs=3
#I want at least 3 datasets (so 3 na's is too many)

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#I'm going to grab just the column with the pvalues and run a multiple comparisons correction using the Benjamini-Hochberg method ("FDR" or "q-value")
tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))   

#unfortunately, the output from that function re-orders the FDR-corrected p-values in order of "significance"
#We would like the FDR-corrected p-values to be in the their original order (i. the order of the rest of our statistical output!). This order is recorded in the index (row numbers) for the p-values:
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])


FalseDiscoveryCorrection(metaOutput)


#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

row.names(metaOutputFDR_OrderbyPval)[c(1:100)]

getwd()

MakeForestPlots("Adm") 
MakeForestPlots("Prl") ##studies with NA omitted from model fitting
MakeForestPlots("Gvin1")
MakeForestPlots("Atp2c2")
MakeForestPlots("Pla2g5")

MakeForestPlots("Tmem125")
MakeForestPlots("Mblac1")
MakeForestPlots("Gm6710")
MakeForestPlots("Btbd19")
MakeForestPlots("Mlycd")
MakeForestPlots("Serpinb1a")
MakeForestPlots("Lrrc3b")
MakeForestPlots("Rdh5")
#Let's see how many results are both statistically significant (using a more conservative FDR) and have a larger estimated Log2FC (>0.10 or <-0.10)

largelogandstatsig <- row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & abs(metaOutputFDR_OrderbyPval[,1])>0.1 & is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
#What if we look at the correlation between datasets again, but only using the top genes found in all 5 datasets?

TopGenesInAll5Datasets<-row.names(metaOutputFDR)[metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>4 & is.na(metaOutputFDR[,6])==FALSE]
##1123


CorMatrixTopGenes<-cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll5Datasets),-1]), use="pairwise.complete.obs")


setwd("C:/Users/annak/Documents/R")
RunMetaAnalysisDiagnostics("Prl")



largelogandstatsig