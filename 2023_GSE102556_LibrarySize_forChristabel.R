
#######################

#GSE102556.1
#This is an RNA-Seq dataset:

rm(SummarizedExperiment)
SummarizedExperiment<-gemma.R::get_dataset_object("GSE102556.1", type = 'se', filter=TRUE, consolidate="average")
SummarizedExperiment
#dim: 24935 78 

#Multi-QC... not available? :(

#Can we get the count data?
https://gemma.msl.ubc.ca/rest/v2/datasets/GSE102556.1/data/raw?quantitationType=Counts
#Yes, counts are available

GSE102556_Counts<-read.delim("19431_GSE102556.1_expmat.unfilt.raw.data.txt", sep="\t", comment.char = "#", header=TRUE, stringsAsFactors = FALSE)
colnames(GSE102556_Counts)

GSE102556_LibrarySize<-apply(GSE102556_Counts[,-c(1:6)], 2, sum)
hist(GSE102556_LibrarySize)
min(GSE102556_LibrarySize)
#[1] 567220.9
#Ouch.
max(GSE102556_LibrarySize)
#[1] 50307736
#That's more respectable

summary(GSE102556_LibrarySize)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 567221 20527044 24493717 24732941 27096051 50307736

50307736/3
#16769245

sum(GSE102556_LibrarySize<16769245)
#[1] 7
#That's not a small number of samples...

GSE102556_LibrarySize

#Does library size relate to brain region at all?
names(GSE102556_LibrarySize)
#Looks like #1-40 are NACC
#And 41-78 are PFC

hist(GSE102556_LibrarySize[c(1:40)])
#Much better
summary(GSE102556_LibrarySize[c(1:40)])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 16039266 20456197 24820923 24623609 28074272 33482446 
#So NACC data is o.k. regarding variation in library size (16-33 million)

GSE102556_LibrarySize[c(1:40)]
