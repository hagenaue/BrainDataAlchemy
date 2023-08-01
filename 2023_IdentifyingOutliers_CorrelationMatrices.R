#Some example code for looking for outlier samples within transcriptional profiling datasets
#Megan Hagenauer, 7-11-2023

#Example: Dataset 81672

#Reading the data in:
GSE81672_Expression<-gemma.R::get_dataset_expression("GSE81672")

#Making a numeric matrix out of the gene expression data columns
GSE81672_Expression_Matrix<-as.matrix(GSE81672_Expression[,-c(1:4)])

#Creating an example sample-sample correlation scatterplot (data for all genes for 1 sample versus the data for all genes for the second sample)
plot(GSE81672_Expression_Matrix[,1]~GSE81672_Expression_Matrix[,2])

#Creating a matrix showing how each sample correlates with every other sample:
GSE81672_CorMatrix<-cor(GSE81672_Expression_Matrix)

#Writing that matrix out to your working directory to save it:
write.csv(GSE81672_CorMatrix, "GSE81672_CorMatrix.csv")

#Creating a hierarchically clustered heatmap illustrating the sample-sample correlation matrix:
heatmap(GSE81672_CorMatrix)

#Creating a boxplot illustrating the sample-sample correlations for each sample. Outliers should be obvious at this point.
boxplot(GSE81672_CorMatrix)
#There is one sample here that looks like a *potential* outlier but it is a little ambiguous



#Here is an Agilent microarray dataset with a more obvious outlier:

GSE84183_Expression<-gemma.R::get_dataset_expression("GSE84183")
GSE84183_Expression_Matrix<-as.matrix(GSE84183_Expression[,-c(1:4)])

#We discovered last year that Agilent datasets often look like they have been log2 transformed twice - this code fixes that for this dataset:

GSE84183_Expression_Matrix_Corrected<-2^GSE84183_Expression_Matrix
hist(GSE84183_Expression_Matrix_Corrected)

#Sample-Sample correlation: scatterplot
plot(GSE84183_Expression_Matrix_Corrected[,1]~GSE84183_Expression_Matrix_Corrected[,2])

GSE84183_CorMatrix<-cor(GSE84183_Expression_Matrix_Corrected)

write.csv(GSE84183_CorMatrix, "GSE84183_CorMatrix.csv")

heatmap(GSE84183_CorMatrix)

boxplot(GSE84183_CorMatrix)
