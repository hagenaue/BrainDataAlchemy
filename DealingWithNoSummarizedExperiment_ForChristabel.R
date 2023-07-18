#Dealing with a lack of Summarized Experiment object in Gemma
#Megan Hagenauer, July 18, 2023

#Dataset: GSE85136 - for Christabel
#This is an RNA-Seq dataset
#It isn't playing well with the summarized experiment or expression set dataset objects

rm(SummarizedExperiment)
SummarizedExperiment<-gemma.R::get_dataset_object("GSE85136", type = 'se')
# Error in `dplyr::mutate()`:
#   ! Can't transform a data frame with `NA` or `""` names.
# Run `rlang::last_trace()` to see where the error occurred.

#But it can be read in just as expression data:
GSE85136_Expression<-gemma.R::get_dataset_expression("GSE85136", filter=TRUE)
str(GSE85136_Expression)
#Classes ‘data.table’ and 'data.frame':	23821 obs. of  28 variables:
colnames(GSE85136_Expression)
#Much of the important metadata is in the column name (sex, stress)
#the colnames are not where the problematic NAs are located

head(row.names(GSE85136_Expression))
#[1] "1" "2" "3" "4" "5" "6"
sum(row.names(GSE85136_Expression)=="")
#[1] 0
#The row names are not where the problematic NAs are located.

GSE85136_Expression[GSE85136_Expression$Probe=="", c(1:4)]
#none
#So the NAs causing problems with reading in the data are not in the Probe column
GSE85136_Expression[GSE85136_Expression$NCBIid=="", c(1:4)]
#quite a few
nrow(GSE85136_Expression[GSE85136_Expression$NCBIid=="", c(1:4)])
#[1] 20
#But not that many.
#Maybe those are the NAs causing issues?

#There isn't an option to consolidate the data representing the same gene for this version of reading in data
#This is RNA-Seq data, so it may not be necessary

max(table(GSE85136_Expression$NCBIid))
#[1] 20
#There do seem to be a few transcripts that are either lack annotation or mapped to the same gene

names(table(GSE85136_Expression$NCBIid)[table(GSE85136_Expression$NCBIid)==20])
#[1] ""
#Looks like these lack NCBI Ids
#Weird. I wonder what the original annotation was for the counts
head(GSE85136_Expression$Probe)
#100009600    100017    100019 100033459 100034251 100034361
head(GSE85136_Expression$NCBIid)
#[1] "100009600" "100017"    "100019"    "100033459" "100034251" "100034361"
#Looks the same to me... huh...

GSE85136_Expression[which(GSE85136_Expression$Probe!=GSE85136_Expression$NCBIid),c(1:4)]
#Looking at the gemma notes it seems like the counts were originally annotated with Ensemble
#When annotated with entrez ids, it looks like sometimes those entrez ids no longer existed or were ambiguous(?)

#So sometimes the NCBIids have two ids and a pipe, e.g.:
#"102633292|102633118"

#Those will be problematic later on
#how many are there?
#grep pulls out the rows of genes including a pattern of string
GSE85136_Expression[grep('\\|', GSE85136_Expression$NCBIid), c(1:4)]
#There are only a handful, we can just toss them out when filtering by row

#Are there other rows that have the same NCBIid?
names(table(GSE85136_Expression$NCBIid)[table(GSE85136_Expression$NCBIid)>1])
#[1] ""          "100039542" "677884" 
#Just a few - they popped up when determining which probes don't match NCBIid as well.


NAsPerColumn<-apply(GSE85136_Expression[,-c(1:4)], 2, function(y) sum(is.na(y)))
max(NAsPerColumn)
#[1] 0
NAsPerRow<-apply(GSE85136_Expression[,-c(1:4)], 1, function(y) sum(is.na(y)))
max(NAsPerRow)
#[1] 0
#No missing gene expression data

#So if we want this to play well with our other code, we'll have to shape it into a SummarizedExperiment object

#Note I accidentally ran some code here to filter the data again and now removed it because it was unnecessary
#Here is just some code renaming the object so the downstream code doesn't break
#Note: Some of the downstream dimensions will be different
GSE85136_Expression_Filtered<-GSE85136_Expression
dim(GSE85136_Expression_Filtered)

#We still need the metadata for this object
colnames(GSE85136_Expression_Filtered)

GSE85136_Samples<-gemma.R::get_dataset_design("GSE85136")
str(GSE85136_Samples)
#'data.frame':	24 obs. of  4 variables:
#Double-checking that it is in the same order as the expression data.frame:
row.names(GSE85136_Samples)
# [1] "ControlMale3"   "ControlFemale2" "GFPFemale1"     "GFPMale1"       "StressFemale2"  "CREMale3"      
# [7] "StressMale1"    "CREFemale1"     "StressFemale1"  "GFPMale3"       "CREFemale2"     "GFPMale2"      
# [13] "StressFemale3"  "StressMale2"    "ControlFemale3" "ControlFemale1" "CREMale2"       "GFPFemale3"    
# [19] "ControlMale2"   "GFPFemale2"     "StressMale3"    "ControlMale1"   "CREMale1"       "CREFemale3"    
colnames(GSE85136_Expression_Filtered)
# [1] "Probe"            "GeneSymbol"       "GeneName"         "NCBIid"           "Control Male 3"   "Control Male 1"  
# [7] "Control Male 2"   "Stress Male 2"    "Stress Male 1"    "Stress Male 3"    "GFP Male 3"       "GFP Male 2"      
# [13] "GFP Male 1"       "CRE Male 3"       "CRE Male 2"       "CRE Male 1"       "Control Female 3" "Control Female 1"
# [19] "Control Female 2" "Stress Female 3"  "Stress Female 2"  "Stress Female 1"  "GFP Female 3"     "GFP Female 2"    
# [25] "GFP Female 1"     "CRE Female 3"     "CRE Female 2"     "CRE Female 1" 

#Interestingly, it is *not* in the same order, and the names don't match.
#To make a version that matches, we would need to remove spaces

ExpressionColnames_NoSpaces<- gsub(" ", "", colnames(GSE85136_Expression_Filtered)[-c(1:4)])
ExpressionColnames_NoSpaces
# [1] "ControlMale3"   "ControlMale1"   "ControlMale2"   "StressMale2"    "StressMale1"    "StressMale3"   
# [7] "GFPMale3"       "GFPMale2"       "GFPMale1"       "CREMale3"       "CREMale2"       "CREMale1"      
# [13] "ControlFemale3" "ControlFemale1" "ControlFemale2" "StressFemale3"  "StressFemale2"  "StressFemale1" 
# [19] "GFPFemale3"     "GFPFemale2"     "GFPFemale1"     "CREFemale3"     "CREFemale2"     "CREFemale1"  

#To get the metadata in the same order, we could use join or merge:
library(plyr)
ExpressionColnames_NoSpaces_toJoin<-data.frame(x=ExpressionColnames_NoSpaces)
str(ExpressionColnames_NoSpaces_toJoin)
GSE85136_Samples_toJoin<-data.frame(x=row.names(GSE85136_Samples), GSE85136_Samples)
str(GSE85136_Samples_toJoin)
GSE85136_Samples_toJoin_inOrder<-join(ExpressionColnames_NoSpaces_toJoin, GSE85136_Samples_toJoin, by="x", type="left")
str(GSE85136_Samples_toJoin_inOrder)
GSE85136_Samples_toJoin_inOrder$x
# [1] "ControlMale3"   "ControlMale1"   "ControlMale2"   "StressMale2"    "StressMale1"    "StressMale3"   
# [7] "GFPMale3"       "GFPMale2"       "GFPMale1"       "CREMale3"       "CREMale2"       "CREMale1"      
# [13] "ControlFemale3" "ControlFemale1" "ControlFemale2" "StressFemale3"  "StressFemale2"  "StressFemale1" 
# [19] "GFPFemale3"     "GFPFemale2"     "GFPFemale1"     "CREFemale3"     "CREFemale2"     "CREFemale1"   
#Now they are in the same order
#I'm not sure if the Summarized Experiment function will be cranky if the names don't match - let's find out.

#Code for making a summarized experiment object
GSE85136_SummarizedExperiment<-SummarizedExperiment(assays=list(counts=GSE85136_Expression_Filtered[,-c(1:4)]), colData=GSE85136_Samples_toJoin_inOrder, rowData=GSE85136_Expression_Filtered[,c(1:4)])
GSE85136_SummarizedExperiment
#dim: 19829 24 
#Not bad...
GSE85136_SummarizedExperiment$genotype
#That works
head(rowData(GSE85136_SummarizedExperiment))
#That seems to work too

GSE85136_ExpressionData_fromSE<-assay(GSE85136_SummarizedExperiment)
str(GSE85136_ExpressionData_fromSE)
#Classes ‘data.table’ and 'data.frame':	19829 obs. of  24 variables:

#Looks good.
#Only issue: I have not consolidated (averaged) by gene.
#That won't really matter much for this dataset, but might for others.
#So we might want to add a coding step for it.


#Taking a snoop at the MultiQC Report:
#https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=9572
#There is *huge* variability in library size, with a whole group of samples having much smaller libraries
#8-74.1 million aligned reads
#All samples GSM2258228-GSM2258239 have library sizes >35 million
#All samples GSM2258240-GSM2258251 have library sizes <13 million 
#Egads
#Downloading the design file from here:
#https://gemma.msl.ubc.ca/experimentalDesign/showExperimentalDesign.html?eeid=9572#
#After sorting by GSM number, it looks like all of the small library size samples are from the wildtype mice
#The samples with the larger library sizes are from the CRE mice (CRE control and CRE Dnmt3a)
#I would toss out all of the CRE mice
#After that, the variation in library size is minimal, although shallow (8-13.7 million)

#Any issue getting counts?
https://gemma.msl.ubc.ca/rest/v2/datasets/GSE85136/data/raw?quantitationType=Counts
#{"error":{"code":404,"message":"Entity with the given identifier does not exist or is not accessible.","errors":{"exceptionMessage":"The identifier was recognised to be 'id', but entity of type 'ubic.gemma.model.common.quantitationtype.QuantitationType' with 'id' equal to 'Counts' does not exist or is not accessible.","exceptionName":"ubic.gemma.web.util.EntityNotFoundException"}},"apiVersion":"2.6.1"}
#Yep - counts not available. Bah.
