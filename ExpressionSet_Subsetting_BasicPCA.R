#Example code for subsetting the dataset and running (basic) principal components analysis
#Megan Hagenauer, July 12 2023

#Gemma has created objects that include both the expression data and the metadata (Expression Set)
#These are compatible with other gene expression analysis code packages

GSE81672_ExpressionSet<-gemma.R::get_dataset_object("GSE81672", type = 'eset')
str(GSE81672_ExpressionSet)

#Example - pulling out the batch information for this dataset:
GSE81672_ExpressionSet[[1]]$block

#... and looking at how many samples fall in each batch:
table(GSE81672_ExpressionSet[[1]]$block)

#Looking at the number of samples present for each of the other variables:
table(GSE81672_ExpressionSet[[1]]$`organism part`)
# basolateral amygdaloid nuclear complex                      nucleus accumbens                      prefrontal cortex 
# 25                                     25                                     25 
# ventral,Ammon's horn 
#                                     24 

table(GSE81672_ExpressionSet[[1]]$phenotype)
# Non-responder,susceptible toward                     resistant to               susceptible toward 
# 28                               16                               12 
# susceptible toward,Responder               wild type genotype 
# 24                               19 

table(GSE81672_ExpressionSet[[1]]$treatment)
# imipramine            ketamine          reference substance role,saline 
# 28                              24                              47 



#Pulling out the numeric expression data for this dataset and putting it in matrix form 
#Similar to what we worked with yesterday - you can run any of the code that we used yesterday with this numeric matrix
library(Biobase)
GSE81672_ExpressionSet_ExpressionDataMatrix<-exprs(GSE81672_ExpressionSet[[1]])
str(GSE81672_ExpressionSet_ExpressionDataMatrix)

#Running basic principal components analysis (PCA)
library(stats)
pca_output<-prcomp(t(GSE81672_ExpressionSet_ExpressionDataMatrix))
#scale.=FALSE as default - after we filter the data, we can make scale=TRUE if we want to emphasize results from genes with less variable expression

PC1<-pca_output$x[,1]
PC2<-pca_output$x[,2]

PC3<-pca_output$x[,3]
PC4<-pca_output$x[,4]

#Output a scree plot for the PCA (no outliers):
png("10 PCA Scree Plot1.png")
plot(summary(pca_output)$importance[2,]~(c(1:length(summary(pca_output)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("10 PCA Scree Plot2.png")
plot(summary(pca_output)$importance[3,]~(c(1:length(summary(pca_output)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis")
dev.off()

#Let's see if these clusters are created due to a variable of interest...
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(GSE81672_ExpressionSet[[1]]$`organism part`))
plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(GSE81672_ExpressionSet[[1]]$`organism part`))


#Subsetting the data:
GSE81672_ExpressionSet_HC<-GSE81672_ExpressionSet[[1]][, GSE81672_ExpressionSet[[1]]$`organism part`=="ventral,Ammon's horn"]
str(GSE81672_ExpressionSet_HC)

table(GSE81672_ExpressionSet_HC$`organism part`)

table(GSE81672_ExpressionSet_HC$treatment)
#                     imipramine                        ketamine reference substance role,saline 
#7                               6                              11 