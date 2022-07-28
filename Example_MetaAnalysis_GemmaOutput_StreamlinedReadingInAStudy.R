#Messing around with meta-analyses of Gemma output
#This file includes the streamlined version of how to read in and extract the preprocessed data and statistical output for a single study (steps #0-5 mentioned below).  
#Megan Hagenauer
#07-25-2022
#Associated R workspace: "~/Documents/Teaching/Grinnell Interns 2022/Workspace_Example_MetaAnalysis_GemmaOutput.RData"


####################

#Goal:

#For practice, I decided to focus on hippocampal chronic stress datasets that weren't included in Yusra's HDRF meta-analysis

#I piggybacked on Jinglin Xiong's efforts triaging chronic stress datasets, using her spreadsheet: "ChronicStressData(after exclusion).xlsx"

#These were the hippocampal datasets that seemed to meet criteria (chronic adult stress):
#GSE59070 - chronic social stress (CSDS), "Ammon's horn", mouse, Agilent platform
#GSE86392 - chronic restraint stress, HC, rat, RNA-Seq
#GSE116009 - CUMS, HC, mouse, RNA-Seq
#GSE109315 -CSDS, ventral (V) HC - mouse, RNA-Seq
#GSE151807 - CMS, HC, mice, Affymetrix 430 2.0

#These datasets were hippocampal, but focused on specific subregions:
#GSE11211 - CSDS, CA1
#GSE56028 - unpredictable chronic mild stress (CUMS), Dentate Gyrus (DG)
#GSE84185 - CUMS, DG
#GSE132819 - CSDS, DG
#GSE84183 - CUMS, DG

#Excluded, included in our HDRF meta-analysis:
#GSE81672 - CSDS, HC - Bagot et al.
#GSE72343 - CSDS, "Ammon's Horn" - Bagot et al.

#Excluded, no associated publication:
#GSE109445 - CUMS, "Ammon's horn"

#Excluded - focused on stress during adolescence instead of adulthood:
#GSE172451- chronic social instability stress, vHC

#Excluded:
#GSE102965 - Learned helplessness, HC - this one may not actually work for the simple type of analysis that we're doing - the Gemma results aren't specific to a single brain area, big frowny face on Gemma

####################

#Overview of general coding steps:

#0) Reading in & visualizing the Log2 Expression data and Metadata for the individual samples - an illustration of where the differential expression results come from. 

#1) Read in the results from Gemma for each dataset

#2) Identify the results for the variables of interest

#3) Remove rows of data that have missing or unambiguous gene annotation

#### Question: Which gene annotation should we use? I was originally planning to use NCBI ID (because it is a little more stable than gene symbol), but if we use gene symbol we can run a half-way decent meta-analysis of data from both rats and mice without having to add in a step where we decode which genes in rats and mice are orthologous using an orthology database, as many mice and rat genes have the same gene symbol (76%, last time I checked).

#4) Collapse results (average) if there is more than row of data representing a single gene 

#5) Extract out the information needed for running a meta-analysis: Use the Log2FC and Tstat to calculate the standard error for the Log2FC, and then use the standard error to calculate the sampling variance.

#6) Combine together the relevant results from different studies into a single data-frame for the effect sizes (Log2FC) and a single data.frame for the sampling variances.  - *This will go in a separate code file.*

#7) Make a correlation matrix to compare the overall results from the different studies. Further visualize the comparison using a hierarchically-clustered heatmap of the correlation matrix. Which studies seem to have similar results? - *This will go in a separate code file.*

#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. - *This will go in a separate code file.*

#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations many times and therefore have a heightened risk of false discovery (false discovery rate correction) - *This will go in a separate code file.*

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. - *This will go in a separate code file.*

#####################################

#Code packages used (may require installation & installation of dependencies):

library(plyr)
library(metafor)


######################################

#0) Reading in & visualizing the Log2 Expression data and Metadata for the individual samples - an illustration of where the differential expression results come from. 

#This section is only slightly stream-lined - it is pretty data-set specific, so not easily functionalized.  
#Therefore, to "streamline" it, I just removed a lot of the extra annotation.

#To play with the Log2 Expression data and Metadata for the individual samples, go to the Gemma website for the dataset, e.g.:
https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=11627

#And then download one of the "profiles" into an intuitively named folder/location:
###the unfiltered version is the Log2 gene expression data for all genes/probes represented on the transcriptional profiling platform for all samples
###the filtered version is the Log2 gene expression data for all genes/probes that meet some quality control criteria (typically sufficient expression levels to be considered reliably "measurable") for all samples

#The file containing the Log2 gene expression data for all samples (in this case, 11627_GSE59070_expmat.data.txt.gz) is in a compressed format and will probably need to be "unzipped" to work with it. Typically double-clicking on the file in your folder will take care of this.

#Then set the working directory for R to the folder where you stashed the downloaded data:
setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/11627_GSE59070_diffExpAnalysis_94802")

#The data is read in as a tab ("\t") delimited text file:
TempDataSet_Filtered<-read.delim("11627_GSE59070_expmat.data.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")
str(TempDataSet_Filtered)

#The first 6 columns are probe annotation - Let's separate them out:
TempDataSet_Filtered_Annotation<-TempDataSet_Filtered[,c(1:6)]
TempDataSet_Log2Expression<-TempDataSet_Filtered[,-c(1:6)]

#The sample meta-data is located in the column names, with individual variables separated by "."
#I would like to separate this data out into a matrix where I can grab the values associated with each of the variables for each sample easily 
#Example of splitting up the column names:
Temp_Colnames_SplitIntoList<-strsplit(colnames(TempDataSet_Log2Expression), "\\.")

#Now each sample (previous column name) is treated as an element in a list, and within that element is all of the sample information. 
#It would be much easier to work with this information if it was in a data.frame format.
TempDataSet_MetaData<-do.call(rbind.data.frame, Temp_Colnames_SplitIntoList)
str(TempDataSet_MetaData)
#The columns don't have variable names yet. 

#From looking at the data.frame, the information that defines our main treatment groups is found in these two columns:
TempDataSet_MetaData[,5]
# [1] acute  acute  acute  acute  acute  acute  8days  8days  8days  8days  8days  8days  13days 13days
# [15] 13days 13days 13days 13days 13days 13days 13days 13days 13days 13days
# Levels: 13days 8days acute
TempDataSet_MetaData[,6]
# [1] Biol        Biol        Biol        Biol        Biol        Biol        Biol        Biol       
# [9] Biol        Biol        Biol        Biol        Biol        Biol        Biol        Biol       
# [17] Biol        Biol        5daysofrest 5daysofrest 5daysofrest 5daysofrest 5daysofrest 5daysofrest
# Levels: 5daysofrest Biol

#I'm going to combine the two columns into a vector containing a single grouping variable:
TempDataSet_StressVar<-paste(TempDataSet_MetaData[,5], TempDataSet_MetaData[,6], sep="_")
str(TempDataSet_StressVar)
#chr [1:24] "acute_Biol" "acute_Biol" "acute_Biol" "acute_Biol" "acute_Biol" "acute_Biol" ...

#Alright, let's practice visualizing the log2 expression data for the samples:
#Let's play with the Log2 gene expression data measured using the first probe listed in the dataset.

#We can visualize the Log2 expression in the different groups using a boxplot:

pdf("Boxplot_Hivep3Expression_vs_StressCondition.pdf", height=5, width=7)
boxplot(as.numeric(TempDataSet_Log2Expression[1,])~TempDataSet_StressVar, ylab="Hivep3 Log2 Expression", xlab="Stress Condition", col=8)
stripchart(as.numeric(TempDataSet_Log2Expression[1,])~TempDataSet_StressVar, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)
dev.off()

#Btw - for future reference, instead of using a row number, you can make charts like this for a particular/gene probe using the gene or probe name too. e.g.,
pdf("Boxplot_Hivep3Expression_vs_StressCondition.pdf", height=5, width=7)
boxplot(as.numeric(TempDataSet_Log2Expression[TempDataSet_Filtered_Annotation$GeneSymbol=="Hivep3",][1,])~TempDataSet_StressVar, ylab="Hivep3 Log2 Expression", xlab="Stress Condition", col=8)
stripchart(as.numeric(TempDataSet_Log2Expression[TempDataSet_Filtered_Annotation$GeneSymbol=="Hivep3",][1,])~TempDataSet_StressVar, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)
dev.off()

#######################################

#1) Read in the differential expression results from Gemma for each dataset:

#To start out with, I downloaded the Gemma differential expression output for the studies of interest.
#I put each study is in its own folder.

#I set working directory to where I downloaded the Gemma differential expression output for this particular study:
setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/11627_GSE59070_diffExpAnalysis_94802")
list.files()

#This procedure for reading in the differential expression results was pretty easily generalizable, so I functionalized it:


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


#Notes about parameters for function: ReadingInGemmaDE()
#ResultSetFileNames should be a character vector containing the names of the result files

#Example Usage of the function:

ReadingInGemmaDE(ResultSetFileNames=c("resultset_ID478782.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	41264 obs. of  13 variables:
# $ Element_Name               : chr  "A_52_P266365" "A_52_P762911" "A_52_P200244" "A_52_P1029978" ...
# $ Gene_Symbol                : chr  "Rbm27" "" "Sptlc2" "Lasp1" ...
# $ Gene_Name                  : chr  "RNA binding motif protein 27" "" "serine palmitoyltransferase, long chain base subunit 2" "LIM and SH3 protein 1" ...
# $ NCBI_ID                    : chr  "225432" "" "20773" "16796" ...
# $ FoldChange_13.d.5.d.of.rest: num  0.03501 -0.1365 -0.000978 0.02177 0.09924 ...
# $ Tstat_13.d.5.d.of.rest     : num  0.371 -0.5066 -0.0117 0.1818 0.2917 ...
# $ PValue_13.d.5.d.of.rest    : num  0.715 0.618 0.991 0.858 0.773 ...
# $ FoldChange_8.d             : num  -0.169 -0.09793 0.06481 -0.000861 -0.1216 ...
# $ Tstat_8.d                  : num  -1.791 -0.3635 0.7737 -0.00719 -0.3576 ...
# $ PValue_8.d                 : num  0.0885 0.7201 0.4482 0.9943 0.7244 ...
# $ FoldChange_13.d            : num  -0.07923 -0.1743 -0.01604 -0.00859 0.04578 ...
# $ Tstat_13.d                 : num  -0.8396 -0.647 -0.1915 -0.0717 0.1346 ...
# $ PValue_13.d                : num  0.411 0.525 0.85 0.944 0.894 ...

####################

#3) Remove rows of data that have missing or unambiguous gene annotation

#Notably, using modern probe annotation, many of these probes now map to more than one gene symbol, and when this happens Gemma separates provides a list of all of the gene symbols that the probe maps to using a pipe. 
#There are also probes that no longer map to any gene symbol (annotated with ""). 

#This data processing step was pretty generalizable, so I functionalized it:

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

#I only want rows annotated with a single Gene Symbol (no pipe):
TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA[-(grep('\\|', TempResultsJoined_NoNA$Gene_Symbol)),]

print("# of rows with good annotation")
print(nrow(TempResultsJoined_NoNA_NoMultimapped))

#For record keeping (sometimes useful for troubleshooting later)
write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_NoNA_NoMultimapped.csv")

rm(TempResultsJoined_NoNA, TempResultsJoined_NoNA_NoMultimapped)

print("Outputted object: TempResultsJoined_NoNA_NoMultimapped")
}

#Example usage of function:

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] 11443
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 11443
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 1838
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 1838
# [1] "# of rows with good annotation"
# [1] 27983
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

str(TempResultsJoined_NoNA_NoMultimapped)
# 'data.frame':	27983 obs. of  13 variables:
#   $ Element_Name               : chr  "A_52_P266365" "A_52_P200244" "A_52_P1029978" "A_52_P674374" ...
# $ Gene_Symbol                : chr  "Rbm27" "Sptlc2" "Lasp1" "Rem1" ...
# $ Gene_Name                  : chr  "RNA binding motif protein 27" "serine palmitoyltransferase, long chain base subunit 2" "LIM and SH3 protein 1" "rad and gem related GTP binding protein 1" ...
# $ NCBI_ID                    : chr  "225432" "20773" "16796" "19700" ...
# $ FoldChange_13.d.5.d.of.rest: num  0.03501 -0.000978 0.02177 -0.000332 -0.161 ...
# $ Tstat_13.d.5.d.of.rest     : num  0.371 -0.01167 0.1818 -0.000529 -0.9336 ...
# $ PValue_13.d.5.d.of.rest    : num  0.715 0.991 0.858 1 0.362 ...
# $ FoldChange_8.d             : num  -0.169 0.06481 -0.000861 0.8238 0.03793 ...
# $ Tstat_8.d                  : num  -1.791 0.7737 -0.00719 1.313 0.22 ...
# $ PValue_8.d                 : num  0.0885 0.4482 0.9943 0.2041 0.8281 ...
# $ FoldChange_13.d            : num  -0.07923 -0.01604 -0.00859 -0.00319 -0.2234 ...
# $ Tstat_13.d                 : num  -0.8396 -0.1915 -0.07173 -0.00509 -1.296 ...
# $ PValue_13.d                : num  0.411 0.85 0.944 0.996 0.21 ...


####################

#2) Identifying the results for the variables of interest for the meta-analysis:

#At this point, we should double-check:

#A) Whether Gemma has used a reference (baseline) group that makes sense for our analysis 

#B) Which group comparisons are actually of interest to us.


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

#Calculating SE:  
TempResultsJoined_NoNA_NoMultimapped_SE<-TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol/TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol

colnames(TempResultsJoined_NoNA_NoMultimapped_SE)<-ComparisonsOfInterest

write.csv(TempResultsJoined_NoNA_NoMultimapped_SE, "TempResultsJoined_NoNA_NoMultimapped_SE.csv")

#For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
#The sampling variance is just the standard error squared.

TempResultsJoined_NoNA_NoMultimapped_SV<-(TempResultsJoined_NoNA_NoMultimapped_SE)^2

colnames(TempResultsJoined_NoNA_NoMultimapped_SV)<-ComparisonsOfInterest

write.csv(TempResultsJoined_NoNA_NoMultimapped_SV, "TempResultsJoined_NoNA_NoMultimapped_SV.csv")

TempMasterResults<-list(Log2FC=TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, Tstat=TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, SE=TempResultsJoined_NoNA_NoMultimapped_SE, SV=TempResultsJoined_NoNA_NoMultimapped_SV)

assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))

print(paste("Output: Named DEResults", GSE_ID, sep="_"))

rm(TempMasterResults, TempResultsJoined_NoNA_NoMultimapped_SV, TempResultsJoined_NoNA_NoMultimapped_SE, TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average, TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_Tstat_Average)

}


#Notes about parameters for function CollapsingDEResults_OneResultPerGene()
#GSE_ID is a string indicating the name of the Gemma dataset
#TempResultsJoined_NoNA_NoMultimapped is the data frame outputted by our previous function
#ComparisonsOfInterest is a character vector containing a list of group comparisons of interest within this dataset. Important: These group comparisons should be listed in exactly the same order as the order that you provide the column names for their associated Fold Change and Tstat output.
#NamesOfFoldChangeColumns is a list containing the columns of TempResultsJoined_NoNA_NoMultimapped containing the FoldChange results for your comparisons of interes, in the same order as the ComparisonsOfInterest vector
#NamesOfTstatColumns is a list containing the columns of TempResultsJoined_NoNA_NoMultimapped containing the Tstat results for your comparisons of interes, in the same order as the ComparisonsOfInterest vector

#Example function usage:

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE59070", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("Stress8days_vs_Acute", "Stress13days_vs_Acute"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_8.d, TempResultsJoined_NoNA_NoMultimapped$FoldChange_13.d), NamesOfTstatColumns=list(TempResultsJoined_NoNA_NoMultimapped$Tstat_8.d, TempResultsJoined_NoNA_NoMultimapped$Tstat_13.d))
# [1] "Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 18503
# [1] "# of rows with unique Gene Symbols:"
# [1] 18503
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 18503     2
# [1] "Output: Named DEResults_GSE59070"

str(DEResults_GSE59070)
# List of 4
# $ Log2FC: num [1:18503, 1:2] -0.00551 0.16795 -0.09224 -0.05196 -0.02993 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18503] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr [1:2] "Stress8days_vs_Acute" "Stress13days_vs_Acute"
# $ Tstat : num [1:18503, 1:2] -0.0319 1.6845 -1.311 -0.3 -0.3194 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18503] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr [1:2] "Stress8days_vs_Acute" "Stress13days_vs_Acute"
# $ SE    : num [1:18503, 1:2] 0.1725 0.0997 0.0704 0.1732 0.0937 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18503] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr [1:2] "Stress8days_vs_Acute" "Stress13days_vs_Acute"
# $ SV    : num [1:18503, 1:2] 0.02977 0.00994 0.00495 0.03 0.00878 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18503] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr [1:2] "Stress8days_vs_Acute" "Stress13days_vs_Acute"

head(DEResults_GSE59070[[1]])

###########

#Alright, let's see how are code handles our other datasets!

#These were the hippocampal datasets that seemed to meet criteria (chronic adult stress):
#GSE59070 - chronic social stress (CSDS), "Ammon's horn", mouse, Agilent platform - Done!
#GSE86392 - chronic restraint stress, HC, rat, RNA-Seq
#GSE116009 - CUMS, HC, mouse, RNA-Seq
#GSE109315 -CSDS, ventral (V) HC - mouse, RNA-Seq
#GSE151807 - CMS, HC, mice, Affymetrix 430 2.0

#Starting with GSE86392:
setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/13263_GSE86392_diffExpAnalysis_105657")

list.files()
# [1] "analysis.results.txt"        "resultset_ID483791.data.txt" "resultset_ID483792.data.txt"
# [4] "resultset_ID483793.data.txt" 

ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID483791.data.txt","resultset_ID483792.data.txt", "resultset_ID483793.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	14001 obs. of  25 variables:
#   $ Element_Name               : chr  "ENSRNOG00000006244" "ENSRNOG00000058658" "ENSRNOG00000003071" "ENSRNOG00000018783" ...
# $ Gene_Symbol                : logi  NA NA NA NA NA NA ...
# $ Gene_Name                  : logi  NA NA NA NA NA NA ...
# $ NCBI_ID                    : logi  NA NA NA NA NA NA ...
# $ FoldChange_fluoxetine      : num  -0.373 -0.412 0.16 -0.313 -0.681 ...
# $ Tstat_fluoxetine           : num  -2.37 -1.06 1.59 -2.95 -2.59 ...
# $ PValue_fluoxetine          : num  0.0556 0.3307 0.1623 0.0257 0.0412 ...
# $ FoldChange_Acupuncture     : num  0.2514 -0.3123 0.1533 -0.1997 -0.0217 ...
# $ Tstat_Acupuncture          : num  1.77 -0.8593 1.554 -1.932 -0.0926 ...
# $ PValue_Acupuncture         : num  0.127 0.423 0.171 0.102 0.929 ...
# $ Gene_Symbol                : logi  NA NA NA NA NA NA ...
# $ Gene_Name                  : logi  NA NA NA NA NA NA ...
# $ NCBI_ID                    : logi  NA NA NA NA NA NA ...
# $ FoldChange_restraint.stress: num  -0.0831 0.2043 -0.0796 0.3783 0.5095 ...
# $ Tstat_restraint.stress     : num  -0.582 0.568 -0.806 3.626 2.027 ...
# $ PValue_restraint.stress    : num  0.5817 0.5908 0.4511 0.011 0.0891 ...
# $ Gene_Symbol                : logi  NA NA NA NA NA NA ...
# $ Gene_Name                  : logi  NA NA NA NA NA NA ...
# $ NCBI_ID                    : logi  NA NA NA NA NA NA ...
# $ FoldChange_pituitary.gland : num  1.118 -0.7316 -0.0253 1.093 0.1882 ...
# $ Tstat_pituitary.gland      : num  8.685 -2.128 -0.3 12.03 0.925 ...
# $ PValue_pituitary.gland     : num  0.000129 0.07744 0.7746 0.00002 0.3906 ...
# $ FoldChange_Ammon.s.horn    : num  -0.00295 -0.1486 -0.1273 -0.08415 -0.8201 ...
# $ Tstat_Ammon.s.horn         : num  -0.0202 -0.4822 -1.487 -0.835 -3.457 ...
# $ PValue_Ammon.s.horn        : num  0.9845 0.6467 0.1875 0.4357 0.0135 ...

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] NA
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] NA
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 0
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

#Hmmmm... that's a problem.
#I'm going to go take a peek at the original result set files and see if there is something screwy here.
#Huh - yeah, oddly enough, that entire dataset lacks either gene symbol or NCBI ID annotation (?)
#There is Ensembl annotation though. I wonder how common that is - we may need to tweak our annotation function.

##########################

#Let's try a different dataset: GSE116009

setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/14942_GSE116009_diffExpAnalysis_118827")

list.files()
[1] "analysis.results.txt"                   "resultset_ID488229.data.txt"           
[3] "resultset_ID488230.data.txt"       

ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID488229.data.txt", "resultset_ID488230.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	18541 obs. of  13 variables:
#   $ Element_Name                                                                               : int  69487 102941 19044 98432 20768 380718 52377 18548 102634056 14829 ...
# $ Gene_Symbol                                                                                : chr  "Ndufaf5" "B630019K06Rik" "Ppox" "Phlpp1" ...
# $ Gene_Name                                                                                  : chr  "NADH:ubiquinone oxidoreductase complex assembly factor 5" "F-box and leucine-rich repeat protein 17 pseudogene" "protoporphyrinogen oxidase" "PH domain and leucine rich repeat protein phosphatase 1" ...
# $ NCBI_ID                                                                                    : chr  "69487" "102941" "19044" "98432" ...
# $ FoldChange_SNCA.human.synuclein.alpha.non.A4.component.of.amyloid.precursor._Overexpression: num  -0.0369 -0.0368 -0.0411 -0.0013 -0.0174 ...
# $ Tstat_SNCA.human.synuclein.alpha.non.A4.component.of.amyloid.precursor._Overexpression     : num  -0.3421 -0.3448 -0.2808 -0.0166 -0.1245 ...
# $ PValue_SNCA.human.synuclein.alpha.non.A4.component.of.amyloid.precursor._Overexpression    : num  0.739 0.737 0.784 0.987 0.903 ...
# $ Gene_Symbol                                                                                : chr  "Ndufaf5" "B630019K06Rik" "Ppox" "Phlpp1" ...
# $ Gene_Name                                                                                  : chr  "NADH:ubiquinone oxidoreductase complex assembly factor 5" "F-box and leucine-rich repeat protein 17 pseudogene" "protoporphyrinogen oxidase" "PH domain and leucine rich repeat protein phosphatase 1" ...
# $ NCBI_ID                                                                                    : chr  "69487" "102941" "19044" "98432" ...
# $ FoldChange_chronic.unpredictable.mild.stress                                               : num  0.0203 0.0643 -0.1131 -0.0916 0.0296 ...
# $ Tstat_chronic.unpredictable.mild.stress                                                    : num  0.19 0.612 -0.783 -1.169 0.214 ...
# $ PValue_chronic.unpredictable.mild.stress                                                   : num  0.853 0.553 0.45 0.267 0.834 ...

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] 0
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 0
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 1
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 1
# [1] "# of rows with good annotation"
# [1] 18540
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

str(TempResultsJoined_NoNA_NoMultimapped)
# 'data.frame':	18540 obs. of  13 variables:
#   $ Element_Name                                                                               : int  69487 102941 19044 98432 20768 380718 52377 18548 102634056 14829 ...
# $ Gene_Symbol                                                                                : chr  "Ndufaf5" "B630019K06Rik" "Ppox" "Phlpp1" ...
# $ Gene_Name                                                                                  : chr  "NADH:ubiquinone oxidoreductase complex assembly factor 5" "F-box and leucine-rich repeat protein 17 pseudogene" "protoporphyrinogen oxidase" "PH domain and leucine rich repeat protein phosphatase 1" ...
# $ NCBI_ID                                                                                    : chr  "69487" "102941" "19044" "98432" ...
# $ FoldChange_SNCA.human.synuclein.alpha.non.A4.component.of.amyloid.precursor._Overexpression: num  -0.0369 -0.0368 -0.0411 -0.0013 -0.0174 ...
# $ Tstat_SNCA.human.synuclein.alpha.non.A4.component.of.amyloid.precursor._Overexpression     : num  -0.3421 -0.3448 -0.2808 -0.0166 -0.1245 ...
# $ PValue_SNCA.human.synuclein.alpha.non.A4.component.of.amyloid.precursor._Overexpression    : num  0.739 0.737 0.784 0.987 0.903 ...
# $ Gene_Symbol                                                                                : chr  "Ndufaf5" "B630019K06Rik" "Ppox" "Phlpp1" ...
# $ Gene_Name                                                                                  : chr  "NADH:ubiquinone oxidoreductase complex assembly factor 5" "F-box and leucine-rich repeat protein 17 pseudogene" "protoporphyrinogen oxidase" "PH domain and leucine rich repeat protein phosphatase 1" ...
# $ NCBI_ID                                                                                    : chr  "69487" "102941" "19044" "98432" ...
# $ FoldChange_chronic.unpredictable.mild.stress                                               : num  0.0203 0.0643 -0.1131 -0.0916 0.0296 ...
# $ Tstat_chronic.unpredictable.mild.stress                                                    : num  0.19 0.612 -0.783 -1.169 0.214 ...
# $ PValue_chronic.unpredictable.mild.stress                                                   : num  0.853 0.553 0.45 0.267 0.834 ...

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE116009", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("CUMS_vs_Control"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_chronic.unpredictable.mild.stress), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_chronic.unpredictable.mild.stress)) 
# [1] "Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 18540
# [1] "# of rows with unique Gene Symbols:"
# [1] 18540
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 18540     1
# [1] "Output: Named DEResults_GSE116009"

str(DEResults_GSE116009)
# List of 4
# $ Log2FC: num [1:18540, 1] 0.08169 0.2029 -0.3167 -0.0491 0.00719 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18540] "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr "CUMS_vs_Control"
# $ Tstat : num [1:18540, 1] 0.5963 0.7636 -0.7863 -0.3708 0.0468 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18540] "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr "CUMS_vs_Control"
# $ SE    : num [1:18540, 1] 0.137 0.266 0.403 0.132 0.153 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18540] "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr "CUMS_vs_Control"
# $ SV    : num [1:18540, 1] 0.0188 0.0706 0.1622 0.0175 0.0235 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18540] "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr "CUMS_vs_Control"

#Looks good!

###################

#GSE109315

setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/18225_GSE109315_diffExpAnalysis_172945")

list.files()
# [3] "analysis.results.txt" 
# [32] "resultset_ID498174.data.txt"                                                                  
# [33] "resultset_ID498175.data.txt" 

ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID498174.data.txt", "resultset_ID498175.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	35191 obs. of  16 variables:
#   $ Element_Name                                          : int  242700 269033 258147 102631609 208760 100039542 102639519 26465 723892 14060 ...
# $ Gene_Symbol                                           : chr  "Ifnlr1" "4930503L19Rik" "Olfr675" "Gm29906" ...
# $ Gene_Name                                             : chr  "interferon lambda receptor 1" "RIKEN cDNA 4930503L19 gene" "olfactory receptor 675" "predicted gene, 29906" ...
# $ NCBI_ID                                               : chr  "242700" "269033" "258147" "102631609" ...
# $ FoldChange_resilient.to.chronic.social.defeat.stress  : num  -2.317 -0.3495 -0.00305 0.7678 -0.425 ...
# $ Tstat_resilient.to.chronic.social.defeat.stress       : num  -2.763 -2.172 -0.667 1.017 -0.606 ...
# $ PValue_resilient.to.chronic.social.defeat.stress      : num  0.0108 0.04 0.511 0.3194 0.5501 ...
# $ FoldChange_susceptible.to.chronic.social.defeat.stress: num  -0.8551 0.05122 -0.00645 0.7547 -0.6661 ...
# $ Tstat_susceptible.to.chronic.social.defeat.stress     : num  -1.129 0.353 -1.566 1.107 -1.053 ...
# $ PValue_susceptible.to.chronic.social.defeat.stress    : num  0.27 0.728 0.13 0.279 0.303 ...
# $ Gene_Symbol                                           : chr  "Ifnlr1" "4930503L19Rik" "Olfr675" "Gm29906" ...
# $ Gene_Name                                             : chr  "interferon lambda receptor 1" "RIKEN cDNA 4930503L19 gene" "olfactory receptor 675" "predicted gene, 29906" ...
# $ NCBI_ID                                               : chr  "242700" "269033" "258147" "102631609" ...
# $ FoldChange_C57BL.6Crl                                 : num  0.9039 1.147 -0.00187 -1.533 -0.6793 ...
# $ Tstat_C57BL.6Crl                                      : num  1.212 8.016 -0.462 -2.284 -1.09 ...
# $ PValue_C57BL.6Crl                                     : num  2.37e-01 3.05e-08 6.48e-01 3.15e-02 2.86e-01 ...

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] 2
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 2
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 8
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 8
# [1] "# of rows with good annotation"
# [1] 35181
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

str(TempResultsJoined_NoNA_NoMultimapped)
# 'data.frame':	35181 obs. of  16 variables:
#   $ Element_Name                                          : int  242700 269033 258147 102631609 208760 100039542 102639519 26465 723892 14060 ...
# $ Gene_Symbol                                           : chr  "Ifnlr1" "4930503L19Rik" "Olfr675" "Gm29906" ...
# $ Gene_Name                                             : chr  "interferon lambda receptor 1" "RIKEN cDNA 4930503L19 gene" "olfactory receptor 675" "predicted gene, 29906" ...
# $ NCBI_ID                                               : chr  "242700" "269033" "258147" "102631609" ...
# $ FoldChange_resilient.to.chronic.social.defeat.stress  : num  -2.317 -0.3495 -0.00305 0.7678 -0.425 ...
# $ Tstat_resilient.to.chronic.social.defeat.stress       : num  -2.763 -2.172 -0.667 1.017 -0.606 ...
# $ PValue_resilient.to.chronic.social.defeat.stress      : num  0.0108 0.04 0.511 0.3194 0.5501 ...
# $ FoldChange_susceptible.to.chronic.social.defeat.stress: num  -0.8551 0.05122 -0.00645 0.7547 -0.6661 ...
# $ Tstat_susceptible.to.chronic.social.defeat.stress     : num  -1.129 0.353 -1.566 1.107 -1.053 ...
# $ PValue_susceptible.to.chronic.social.defeat.stress    : num  0.27 0.728 0.13 0.279 0.303 ...
# $ Gene_Symbol                                           : chr  "Ifnlr1" "4930503L19Rik" "Olfr675" "Gm29906" ...
# $ Gene_Name                                             : chr  "interferon lambda receptor 1" "RIKEN cDNA 4930503L19 gene" "olfactory receptor 675" "predicted gene, 29906" ...
# $ NCBI_ID                                               : chr  "242700" "269033" "258147" "102631609" ...
# $ FoldChange_C57BL.6Crl                                 : num  0.9039 1.147 -0.00187 -1.533 -0.6793 ...
# $ Tstat_C57BL.6Crl                                      : num  1.212 8.016 -0.462 -2.284 -1.09 ...
# $ PValue_C57BL.6Crl                                     : num  2.37e-01 3.05e-08 6.48e-01 3.15e-02 2.86e-01 ...

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE109315", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest = c("StressResilient_Vs_Control", "StressSusceptible_Vs_Control"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_resilient.to.chronic.social.defeat.stress, TempResultsJoined_NoNA_NoMultimapped$FoldChange_susceptible.to.chronic.social.defeat.stress), NamesOfTstatColumns=list(TempResultsJoined_NoNA_NoMultimapped$Tstat_resilient.to.chronic.social.defeat.stress, TempResultsJoined_NoNA_NoMultimapped$Tstat_susceptible.to.chronic.social.defeat.stress))
# [1] "Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 35181
# [1] "# of rows with unique Gene Symbols:"
# [1] 35180
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 35180     2
# [1] "Output: Named DEResults_GSE109315"

str(DEResults_GSE109315)
# List of 4
# $ Log2FC: num [1:35180, 1:2] -0.789 0.221 0.392 -0.725 -0.472 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:35180] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# .. ..$ : chr [1:2] "StressResilient_Vs_Control" "StressSusceptible_Vs_Control"
# $ Tstat : num [1:35180, 1:2] -0.903 0.575 2.195 -1.009 -2.078 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:35180] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# .. ..$ : chr [1:2] "StressResilient_Vs_Control" "StressSusceptible_Vs_Control"
# $ SE    : num [1:35180, 1:2] 0.874 0.384 0.179 0.718 0.227 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:35180] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# .. ..$ : chr [1:2] "StressResilient_Vs_Control" "StressSusceptible_Vs_Control"
# $ SV    : num [1:35180, 1:2] 0.7637 0.1475 0.0319 0.516 0.0515 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:35180] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# .. ..$ : chr [1:2] "StressResilient_Vs_Control" "StressSusceptible_Vs_Control"

############################

#GSE151807

setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/17922_GSE151807_diffExpAnalysis_172819")
list.files()
# [1] "analysis.results.txt" 
# [7] "resultset_ID498007.data.txt"  

ReadingInGemmaDE(ResultSetFileNames=c("resultset_ID498007.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	45089 obs. of  7 variables:
# $ Element_Name                  : chr  "1460444_at" "1447047_at" "1460302_at" "1434792_at" ...
# $ Gene_Symbol                   : chr  "Arrb1" "" "Thbs1" "2010320M18Rik" ...
# $ Gene_Name                     : chr  "arrestin, beta 1" "" "thrombospondin 1" "RIKEN cDNA 2010320M18 gene" ...
# $ NCBI_ID                       : chr  "109689" "" "21825" "72093" ...
# $ FoldChange_chronic.mild.stress: num  -0.00673 0.979 -0.07901 -0.1962 0.1982 ...
# $ Tstat_chronic.mild.stress     : num  -0.276 8.843 -0.598 -2.053 3.306 ...
# $ PValue_chronic.mild.stress    : num  0.796 0.000903 0.5819 0.1093 0.02978 ...

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] 13885
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 13885
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 659
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 659
# [1] "# of rows with good annotation"
# [1] 30545
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

str(TempResultsJoined_NoNA_NoMultimapped)
# 'data.frame':	30545 obs. of  7 variables:
#   $ Element_Name                  : chr  "1460444_at" "1460302_at" "1434792_at" "1424623_at" ...
# $ Gene_Symbol                   : chr  "Arrb1" "Thbs1" "2010320M18Rik" "Serpinb5" ...
# $ Gene_Name                     : chr  "arrestin, beta 1" "thrombospondin 1" "RIKEN cDNA 2010320M18 gene" "serine (or cysteine) peptidase inhibitor, clade B, member 5" ...
# $ NCBI_ID                       : chr  "109689" "21825" "72093" "20724" ...
# $ FoldChange_chronic.mild.stress: num  -0.00673 -0.07901 -0.1962 0.04144 -0.08849 ...
# $ Tstat_chronic.mild.stress     : num  -0.276 -0.598 -2.053 1.121 -1.004 ...
# $ PValue_chronic.mild.stress    : num  0.796 0.582 0.109 0.325 0.372 ...

#Note: I came back to this dataset later because I was concerned that it had read in incorrectly because so many of the results had bizarrely small SE associated with them.
#It turns out that the dataset has a many results with suspiciously large t-stats - esp. for a dataset with n=3/group!
hist(TempResultsJoined_NoNA_NoMultimapped$Tstat_chronic.mild.stress, breaks=100)

summary(TempResultsJoined_NoNA_NoMultimapped$Tstat_chronic.mild.stress)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -41.980  -1.730  -0.293  -0.258   1.162  71.490 
#Extremely suspicious

summary(TempResultsJoined_NoNA_NoMultimapped$PValue_chronic.mild.stress)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000002 0.0476000 0.2173000 0.3164795 0.5399000 1.0000000 
#Over 1/4 of the genes have results with p<0.05

#To double-check that this wasn't a problem with my processing of the data, I went back and double-checked the resultset in Excel, and there is indeed a large number of extreme T-stats. I also re-downloaded the data to double-check it, and it looked the same. 


CollapsingDEResults_OneResultPerGene(GSE_ID="GSE151807", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest = c("CMS_vs_Ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_chronic.mild.stress), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_chronic.mild.stress))
# [1] "Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 18337
# [1] "# of rows with unique Gene Symbols:"
# [1] 18337
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 18337     1
# [1] "Output: Named DEResults_GSE151807"

str(DEResults_GSE151807)
# List of 4
# $ Log2FC: num [1:18337, 1] -0.2051 -0.1937 -0.2101 -0.0245 -0.3498 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18337] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr "CMS_vs_Ctrl"
# $ Tstat : num [1:18337, 1] -3.18 -11.76 -2.5 -0.91 -5.46 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18337] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr "CMS_vs_Ctrl"
# $ SE    : num [1:18337, 1] 0.0645 0.0165 0.0841 0.0269 0.0641 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18337] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr "CMS_vs_Ctrl"
# $ SV    : num [1:18337, 1] 0.004157 0.000271 0.007074 0.000722 0.00411 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:18337] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010F05Rik" ...
# .. ..$ : chr "CMS_vs_Ctrl"



###################

#Alright!  On to cross-dataset analysis!  I'll put that code in a new file...

#Before doing anything more though, I am going to save my workspace:

save.image("~/Documents/Teaching/Grinnell Interns 2022/Workspace_Example_MetaAnalysis_GemmaOutput.RData")

#And take some quick notes on the code packages that are loaded & their versions:
#(this information is often needed for publication, and can also be helpful for debugging if someone else tries to reuse my code and it doesn't work):

sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] metafor_2.0-0    Matrix_1.2-10    compute.es_0.2-5 plyr_1.8.4      
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.4.6     lattice_0.20-35  crayon_1.3.4     dplyr_0.8.0.1    assertthat_0.2.1
# [6] grid_3.4.1       R6_2.4.0         nlme_3.1-131     magrittr_1.5     pillar_1.3.1    
# [11] rlang_0.3.1      rstudioapi_0.6   tools_3.4.1      glue_1.3.1       purrr_0.3.2     
# [16] compiler_3.4.1   pkgconfig_2.0.2  tidyselect_0.2.5 tibble_2.1.1  


##############

#The CMS dataset with the suspiciously large t-stats (despite super small sample size) is really throwing off the meta-analysis.
#I'm curious to know how many more datasets would be required to dilute it's effect on the results, so I'm going to toss in the HDRF datasets too.

#Included in our HDRF meta-analysis:
#GSE81672 - CSDS, HC - Bagot et al.
#GSE72343 - CSDS, "Ammon's Horn" - Bagot et al.


#GSE81672

setwd("~/Documents/SideProjects/BrainGMT/Gemma/Hippocampus/16574_GSE81672_diffExpAnalysis_158324")
list.files()
# [1] "analysis.results.txt"  
# [23] "resultset_ID492708.data.txt"                                                         
# [24] "resultset_ID492709.data.txt"  

ReadingInGemmaDE(ResultSetFileNames=c("resultset_ID492708.data.txt", "resultset_ID492709.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	21770 obs. of  25 variables:
#   $ Element_Name                               : int  71673 13383 102638195 110253 22334 26383 226090 12877 21899 66158 ...
# $ Gene_Symbol                                : chr  "Rnf215" "Dlg1" "2700003A03Rik" "Triobp" ...
# $ Gene_Name                                  : chr  "ring finger protein 215" "discs large MAGUK scaffold protein 1" "RIKEN cDNA 2700003A03 gene" "TRIO and F-actin binding protein" ...
# $ NCBI_ID                                    : chr  "71673" "13383" "102638195" "110253" ...
# $ FoldChange_susceptible.toward_Responder    : num  -0.02669 -0.00939 -0.6363 0.1519 0.03877 ...
# $ Tstat_susceptible.toward_Responder         : num  -0.438 -0.191 -1.163 1.486 0.877 ...
# $ PValue_susceptible.toward_Responder        : num  0.667 0.851 0.26 0.155 0.392 ...
# $ FoldChange_resistant.to                    : num  -0.0424 0.0136 -1.047 0.0211 0.0529 ...
# $ Tstat_resistant.to                         : num  -0.549 0.219 -1.509 0.163 0.944 ...
# $ PValue_resistant.to                        : num  0.59 0.829 0.149 0.873 0.358 ...
# $ FoldChange_susceptible.toward              : num  -0.0736 -0.0214 0.1094 0.0223 -0.0696 ...
# $ Tstat_susceptible.toward                   : num  -0.882 -0.318 0.146 0.159 -1.15 ...
# $ PValue_susceptible.toward                  : num  0.39 0.754 0.886 0.875 0.265 ...
# $ FoldChange_Non.responder_susceptible.toward: logi  NA NA NA NA NA NA ...
# $ Tstat_Non.responder_susceptible.toward     : logi  NA NA NA NA NA NA ...
# $ PValue_Non.responder_susceptible.toward    : logi  NA NA NA NA NA NA ...
# $ Gene_Symbol                                : chr  "Rnf215" "Dlg1" "2700003A03Rik" "Triobp" ...
# $ Gene_Name                                  : chr  "ring finger protein 215" "discs large MAGUK scaffold protein 1" "RIKEN cDNA 2700003A03 gene" "TRIO and F-actin binding protein" ...
# $ NCBI_ID                                    : chr  "71673" "13383" "102638195" "110253" ...
# $ FoldChange_imipramine                      : num  -0.00816 -0.05569 0.04738 -0.0571 -0.05889 ...
# $ Tstat_imipramine                           : num  -0.111 -0.941 0.072 -0.465 -1.108 ...
# $ PValue_imipramine                          : num  0.912 0.359 0.943 0.648 0.282 ...
# $ FoldChange_ketamine                        : num  -0.02506 0.00422 -0.00965 0.1859 -0.02469 ...
# $ Tstat_ketamine                             : num  -0.3262 0.0681 -0.014 1.442 -0.4433 ...
# $ PValue_ketamine                            : num  0.748 0.947 0.989 0.166 0.663 ...

FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] "# of rows with missing NCBI annotation:"
# [1] 0
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 0
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 4
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 4
# [1] "# of rows with good annotation"
# [1] 21766
# [1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

str(TempResultsJoined_NoNA_NoMultimapped)
#'data.frame':	21766 obs. of  25 variables:

CollapsingDEResults_OneResultPerGene(GSE_ID="GSE81672", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest = c("StressResistent_vs_Ctrl", "StressSusceptible_vs_Ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_resistant.to, TempResultsJoined_NoNA_NoMultimapped$FoldChange_susceptible.toward), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_resistant.to, TempResultsJoined_NoNA_NoMultimapped$Tstat_susceptible.toward ))
# [1] "Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:"
# [1] "# of rows with unique NCBI IDs:"
# [1] 21766
# [1] "# of rows with unique Gene Symbols:"
# [1] 21766
# [1] "Dimensions of Fold Change matrix, averaged by gene symbol:"
# [1] 21766     2
# [1] "Output: Named DEResults_GSE81672"

str(DEResults_GSE81672)
# List of 4
# $ Log2FC: num [1:21766, 1:2] 0.155 0.1 0.7492 -0.1075 -0.0212 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:21766] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# .. ..$ : chr [1:2] "StressResistent_vs_Ctrl" "StressSusceptible_vs_Ctrl"
# $ Tstat : num [1:21766, 1:2] 0.21 1.025 1.499 -0.984 -0.258 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:21766] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# .. ..$ : chr [1:2] "StressResistent_vs_Ctrl" "StressSusceptible_vs_Ctrl"
# $ SE    : num [1:21766, 1:2] 0.7395 0.0975 0.4998 0.1093 0.082 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:21766] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# .. ..$ : chr [1:2] "StressResistent_vs_Ctrl" "StressSusceptible_vs_Ctrl"
# $ SV    : num [1:21766, 1:2] 0.54687 0.00951 0.2498 0.01194 0.00673 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:21766] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# .. ..$ : chr [1:2] "StressResistent_vs_Ctrl" "StressSusceptible_vs_Ctrl"

####################

#GSE72343

setwd("~/Documents/SideProjects/BrainGMT/Gemma/NucleusAccumbens/16860_GSE72343_diffExpAnalysis_159213")
list.files()
# [1] "analysis.results.txt"  
# [89] "resultset_ID493087.data.txt"                                                                                            
# [90] "resultset_ID493088.data.txt"                                                                                            
# [91] "resultset_ID493089.data.txt"  

ReadingInGemmaDE(ResultSetFileNames=c("resultset_ID493087.data.txt" , "resultset_ID493088.data.txt", "resultset_ID493089.data.txt"))
#[1] "Outputted object: TempResultsJoined"

str(TempResultsJoined)
# 'data.frame':	22948 obs. of  31 variables:
#   $ Element_Name                                                          : int  26440 330721 72183 14027 67685 218581 13483 105651 242608 28088 ...
# $ Gene_Symbol                                                           : chr  "Psma1" "Nek5" "Snx6" "Evpl" ...
# $ Gene_Name                                                             : chr  "proteasome (prosome, macropain) subunit, alpha type 1" "NIMA (never in mitosis gene a)-related expressed kinase 5" "sorting nexin 6" "envoplakin" ...
# $ NCBI_ID                                                               : chr  "26440" "330721" "72183" "14027" ...
# $ FoldChange_amygdala                                                   : num  0.118 -2.008 -0.0371 0.8377 -0.2657 ...
# $ Tstat_amygdala                                                        : num  3.7 -8.9 -1.67 13.87 -5.58 ...
# $ PValue_amygdala                                                       : num  3.26e-04 5.55e-15 9.70e-02 0.00 1.43e-07 ...
# $ FoldChange_nucleus.accumbens                                          : num  0.09255 -1.205 -0.1211 1.103 0.00834 ...
# $ Tstat_nucleus.accumbens                                               : num  2.898 -5.344 -5.453 18.25 0.175 ...
# $ PValue_nucleus.accumbens                                              : num  4.44e-03 4.20e-07 2.58e-07 0.00 8.61e-01 ...
# $ FoldChange_Ammon.s.horn                                               : num  0.173 -0.566 0.015 -0.998 0.644 ...
# $ Tstat_Ammon.s.horn                                                    : num  5.425 -2.51 0.678 -16.52 13.52 ...
# $ PValue_Ammon.s.horn                                                   : num  2.92e-07 1.34e-02 4.99e-01 0.00 0.00 ...
# $ Gene_Symbol                                                           : chr  "Psma1" "Nek5" "Snx6" "Evpl" ...
# $ Gene_Name                                                             : chr  "proteasome (prosome, macropain) subunit, alpha type 1" "NIMA (never in mitosis gene a)-related expressed kinase 5" "sorting nexin 6" "envoplakin" ...
# $ NCBI_ID                                                               : chr  "26440" "330721" "72183" "14027" ...
# $ FoldChange_28.d_chronic.social.defeat.stress_1.h_aggressor.re.exposure: num  -0.1401 -0.00876 -0.09446 0.2622 -0.1254 ...
# $ Tstat_28.d_chronic.social.defeat.stress_1.h_aggressor.re.exposure     : num  -4.901 -0.0434 -4.749 4.847 -2.942 ...
# $ PValue_28.d_chronic.social.defeat.stress_1.h_aggressor.re.exposure    : num  2.92e-06 9.66e-01 5.55e-06 3.67e-06 3.89e-03 ...
# $ FoldChange_28.d_chronic.social.defeat.stress                          : num  -0.1368 0.0866 -0.0527 0.1721 -0.2152 ...
# $ Tstat_28.d_chronic.social.defeat.stress                               : num  -4.784 0.429 -2.648 3.181 -5.048 ...
# $ PValue_28.d_chronic.social.defeat.stress                              : num  4.79e-06 6.69e-01 9.14e-03 1.86e-03 1.55e-06 ...
# $ Gene_Symbol                                                           : chr  "Psma1" "Nek5" "Snx6" "Evpl" ...
# $ Gene_Name                                                             : chr  "proteasome (prosome, macropain) subunit, alpha type 1" "NIMA (never in mitosis gene a)-related expressed kinase 5" "sorting nexin 6" "envoplakin" ...
# $ NCBI_ID                                                               : chr  "26440" "330721" "72183" "14027" ...
# $ FoldChange_resilient                                                  : num  -0.0156 -0.2394 -0.0251 0.0333 -0.0838 ...
# $ Tstat_resilient                                                       : num  -0.564 -1.226 -1.305 0.636 -2.032 ...
# $ PValue_resilient                                                      : num  0.5736 0.2226 0.1942 0.526 0.0442 ...
# $ FoldChange_susceptible                                                : num  0.0122 -0.1789 0.012 -0.0282 -0.0128 ...
# $ Tstat_susceptible                                                     : num  0.441 -0.916 0.625 -0.539 -0.311 ...
# $ PValue_susceptible                                                    : num  0.66 0.361 0.533 0.591 0.756 ...

#Hmmm.... It turns out that this analysis isn't broken down by brain region.
#If we wanted to use this dataset we would need to re-analyze it using only the hippocampal data.

