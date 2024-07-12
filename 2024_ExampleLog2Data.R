#Quickly checking on the range of the gene expression data for each dataset
#Megan Hagenauer
#2024-07-12

#We can download the processed expression data for any particular dataset using this code:
Expression<-gemma.R::get_dataset_processed_expression("GSE81672")
str(Expression)
#Classes ‘data.table’ and 'data.frame':	35205 obs. of  103 variables:
#The first four columns are row metadata: Probe, GeneSymbol, GeneName, NCBIid
#The rest of the columns are gene expression values for each subject
#For reference, this is an RNA-Seq dataset

#You can visualize the distribution of the gene expression data for any particular sample using a histogram:
hist(Expression$`Sample 32: NAC_susceptible_saline`)
#You can see that the range of the x-axis is around -5 to 12
#These are log2 gene expression values - in this case, log2 counts per million (cpm)
#The large spike on the left side of the histogram ("floor effect") are all of the genes that aren't truly expressed or have too low of expression to be measurable

#If you want to visualize the distribution of gene expression data for the entire study, you will need to grab all of the columns that aren't row metadata (i.e., not rows 1-4) and force them into the format of a numeric matrix first:
hist(as.matrix(Expression[,-c(1:4)]))

#We can make this pretty by adding a title and x-axis label:
hist(as.matrix(Expression[,-c(1:4)]), main="Histogram", xlab="Log2 Expression")

#You can also change the color and scaling:
hist(as.matrix(Expression[,-c(1:4)]), main="Histogram", xlab="Log2 Expression", col="red", cex.axis=1.3, cex.lab=1.3)

#You can save the histogram using "export" in the Plots window
#We can also automatically save the histogram by outputting it as a graphics file:
pdf("HistogramForGSE81672.pdf", height=4, width=4)
hist(as.matrix(Expression[,-c(1:4)]), main="Histogram", xlab="Log2 Expression", col="red")
dev.off()
#That will write out into your working directory
#If you don't know where that is, you can find out using:
getwd()

#We can also pull out numeric values summarizing the distribution, e.g.:
min(as.matrix(Expression[,-c(1:4)]))
#[1] -5.8601
median(as.matrix(Expression[,-c(1:4)]))
#[1] -2.1651
max(as.matrix(Expression[,-c(1:4)]))
#[1] 12.312
#Or to get more of an overview:
summary(as.matrix(Expression[,-c(1:4)]))

#If we make boxplots for each subject, we see that the data has been normalized:
boxplot(Expression[,-c(1:4)])
#all of the boxes have the same range


#...so why are we doing this?

#Here's an Agilent dataset that we've found to have issues before...

Expression<-gemma.R::get_dataset_processed_expression("GSE84183")
str(Expression)
#Classes ‘data.table’ and 'data.frame':	59305 obs. of  68 variables:
hist(as.matrix(Expression[,-c(1:4)]))
min(as.matrix(Expression[,-c(1:4)]))
#[1] 2.2011
median(as.matrix(Expression[,-c(1:4)]))
#[1] 2.2975
max(as.matrix(Expression[,-c(1:4)]))
#[1] 4.2235
#The range is between 2 and 4.5
#Why is that?

#The range previously was between -5 to 12 - that was a log2 RNA-Seq dataset
#Log2 Microarray datasets often have a range between 4-15
log2(4)
#[1] 2
log2(8)
#[1] 3
log2(12)
#[1] 3.584963
log2(15)
#[1] 3.906891
#This dataset has been log2 transformed twice!

#If we reverse the log2 transformation by performing 2^X we get a normal-looking distribution:
hist(2^as.matrix(Expression[,-c(1:4)]))


#We need to double check all of our datasets for this problem - especially the datasets lacking raw data, with Agilent microarray data being the main culprit

#We could do this one at a time
#... or we could make a loop to do the operation for *all* of our datasets


