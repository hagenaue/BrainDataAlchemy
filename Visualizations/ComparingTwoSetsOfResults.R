#Example Code for Comparing Two Sets of Differential Expression Results
#Megan Hagenauer
#2025-06-11

###############

#Reading in the differential expression results for the two datasets:

#Set your working directory to where your differential expression results are stored
#This can be done using the GUI: "Session"->"Set Working Directory"


#To read in the results, change the file names below to your actual file names:
Dataset1_DE<-read.csv("Dataset1_DEResults.csv", header=TRUE, stringsAsFactors = FALSE)

str(Dataset1_DE)

Dataset2_DE<-read.csv("Dataset2_DEResults.csv", header=TRUE, stringsAsFactors = FALSE)

###############

#Preparing to join/merge datasets:

#Before joining the results, we will need to make sure that the columns representing the results from each dataset are distinctly named, e.g.,

#Renaming the columns for Dataset 1:

#First, peeking at the original column names:
colnames(Dataset1_DE)

colnames(Dataset1_DE)<-paste("Dataset1", colnames(Dataset1_DE), sep="_")

#Double check that they are renamed:
colnames(Dataset1_DE)

#Renaming the columns for Dataset 2:

#First, peeking at the original column names:
colnames(Dataset2_DE)

colnames(Dataset2_DE)<-paste("Dataset2", colnames(Dataset2_DE), sep="_")

#Double check that they are renamed:
colnames(Dataset2_DE)


#To compare the results, we will need to align them by some sort of gene identifier, so that we can create a data frame with each row representing a gene, and the columns representing the results from each of the studies.

#To do this, we will need to figure out which column in each dataset contains the same type of gene identifier (e.g. Gene Symbol or Entrez ID) as what is found in the other dataset:

colnames(Dataset1_DE)
colnames(Dataset2_DE)

#In this example, both datasets have Gene Symbol annotation, and that annotation is found in column1 in one dataset and column13 in the other
#Note: this process is more complicated if you are joining results from two different species - in that case, we will need to identify homologous gene identifiers in the two datasets using a homolog/ortholog database


#The gene identifier column will need to be renamed so that the name is identical in both datasets, e.g.,

colnames(Dataset1_DE)[1]<-"GeneSymbol"
colnames(Dataset2_DE)[13]<-"GeneSymbol"

#Also important: If there are duplicated gene identifiers in either dataset, you will get the combination of both when joining
#e.g., if DataSet1 has two rows labeled Bdnf and Dataset2 has three rows labeled Bdnf, we will end up with 6 rows in the joined dataset representing each combination of the Bdnf rows:
# Bdnf1-Bdnf1
# Bdnf1-Bdnf2
# Bdnf1-Bdnf3
# Bdnf2-Bdnf1
# Bdnf2-Bdnf2
# Bdnf2-Bdnf3

#... so ideally there should either not be many duplicated gene identifiers, or only have duplicated gene identifiers in one dataset so that you don't have a combinatorial explosion...

#Let's check for that:
sum(duplicated(Dataset1_DE$GeneSymbol))
#[1] 18
sum(duplicated(Dataset2_DE$GeneSymbol))
#[1] 0
#Not bad. I can run with that.

############

#Joining datasets:

install.packages("plyr")
library(plyr)

#To join the two datasets, we will need to tell R which column contains the gene identifier in both datasets that we are using for the alignment
#In this example, that column is named "GeneSymbol"
Dataset1_vs_Dataset2<-join(Dataset1_DE,Dataset2_DE, by="GeneSymbol", type="inner")
#I chose the type "inner" so that my output only contains rows of genes that are found in both datasets
#Other options are "left" (all genes in Dataset1, regardless of whether they are in Dataset2), "right" (the opposite), and "full" (all genes that are in either dataset)
#If you choose those options, there will just be NAs whenever a gene is not found in the results for those datasets

#Side note:
#Over the years, I've seen some folks have trouble with the join function
#The merge function works similarly and can be substituted

#Examining the joined dataframe:
str(Dataset1_vs_Dataset2)

#I would write this out:
write.csv(Dataset1_vs_Dataset2, "Dataset1_vs_Dataset2.csv")

#If you want to make a pretty table showing the results from both datasets for the top genes as identified in one Dataset (e.g., FDR<0.05), that may be easier to do in Excel using the sort function and formatting functions
#... but you can also code it in R if you're feeling fancy...

##########################

#Visualizing comparisons: Scatterplots

#We can make scatterplots comparing results

#e.g., Comparing the Log2FC for all genes:
#Replace "Dataset1_Log2FC" and "Dataset2_Log2FC" with the correct column names

pdf("Scatterplot_Dataset1_vs_Dataset2_Log2FC_AllGenes.pdf", height=5, width=5)
plot(Dataset2_Log2FC~Dataset1_Log2FC, data=Dataset1_vs_Dataset2, xlab="Dataset1 Log2FC", ylab="Dataset2 Log2FC")
TrendLine<-lm(Dataset2_Log2FC~Dataset1_Log2FC, data=Dataset1_vs_Dataset2) #Fitting a trendline
abline(TrendLine, col="red", lwd=3)
dev.off()

#This will give you the linear regression output for that relationship:
summary.lm(TrendLine)
rm(TrendLine)

#If you are using all of the genes represented in the two datasets, you will probably find that the R2 is very small (<0.10) but the p-value is really strong (because there will typically be 12,000-16,000 datapoints...)

#Calculating the Pearson rank correlation:
cor.test(Dataset1_vs_Dataset2$Dataset2_Log2FC, Dataset1_vs_Dataset2$Dataset1_Log2FC, method="pearson")

#Calculating the Spearman rank correlation:
cor.test(Dataset1_vs_Dataset2$Dataset2_Log2FC, Dataset1_vs_Dataset2$Dataset1_Log2FC, method="spearman")

##If you want to make a scatterplot or calculate the linear regression relationship between the Log2FC for just the significant genes from one dataset, the only thing that would need to change in the code above is the information about the data source in the plot() and lm() functions
# The data source would instead be the dataset subsetted down to the significant genes
#For example, this would subset the data down to the rows of genes with an FDR<0.05 in Dataset1:
data=Dataset1_vs_Dataset2[Dataset1_vs_Dataset2$Dataset1_FDR<0.05,]
#to make the code work, you would need to replace that with the correct column name

##########################

#Sometimes folks make scatterplots with the t-statistics
#Especially if the data are noisy
#And RRHOs often use T-statistics
#...but Tstatistics are often missing in differential expression output

#Here's how to calculate them:
#To make this code work, replace with column names containing the Log2FC and Standard Error (SE) for each dataset

Dataset1_vs_Dataset2$Dataset1_Tstat<-(Dataset1_vs_Dataset2$Dataset1_Log2FC/Dataset1_vs_Dataset2$Dataset1_SE)

Dataset1_vs_Dataset2$Dataset2_Tstat<-(Dataset1_vs_Dataset2$Dataset2_Log2FC/Dataset1_vs_Dataset2$Dataset2_SE)

##########################

#Visualizing comparisons: RRHO

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("RRHO")

library(RRHO)
#For making individual RRHO plots/stats:

TempDF<-Dataset1DE_vs_Dataset2

dim(TempDF)

#Remove rows that lack stats output in either dataset:
#e.g.,

TempDF<-TempDF[is.na(TempDF$Dataset1_Log2FC)==FALSE & is.na(TempDF$Dataset2_Log2FC)==FALSE,]

#Determining how many genes are present in the dataframe once the rows lacking stats are removed:
dim(TempDF)

#There also needs to be a unique gene identifier for each row
#Sometimes this causes errors if there are duplicated gene identifiers 
#RRHO asks for these identifiers because I believe it doesn't assume that the input for the two lists is already aligned so it runs its own join or merge function
#But since our datasets are already joined, this can actually just be the row index:
TempDF$RowIndex<-c(1:nrow(TempDF))

list1<-data.frame(GeneRow=TempDF$RowIndex, Metric=TempDF$Dataset1_Tstat)

list2<-data.frame(GeneRow=TempDF$RowIndex, Metric=TempDF$Dataset2_Tstat)

#You'll want to change this code to your desired output directory:
RRHO(list1, list2, labels=c("Dataset1 T-stat", "Dataset2 T-stat"), plots=TRUE, alternative="two.sided", outputdir="~/Dataset2_Archs4/RRHO_AllGenes", BY=TRUE, log10.ind=TRUE)

#Note: the colors in the RRHO and p-values invert for a negative correlation. 
#A little awkward.

#Getting rank-based stats to accompany the RRHO:
#spearman rank correlation:

cor.test(list1[,2], list2[,2], method="spearman")
