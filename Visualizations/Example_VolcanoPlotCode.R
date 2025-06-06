#Example Volcano Plot Code
#Megan Hagenauer
#June 6 2025

################

#Set the working directory to where you are storing the R output files for your project
#This can be done using coding or using the drop-down menu in R studio ("Session" ->"Set Working Directory").
#If you use the drop down menu, I recommend grabbing the code that appears in the console and saving it in your code so that you have documentation as to where your input file came from.

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2024_SophieMensch_ADHD/ROutput_And_Results")

#############

#Note: start with an empty workspace - some of these objects are named similarly to objects that we've used in other coding tasks

rm(list = ls())

#############

#Read in meta-analysis results (or any other differential expression results):

Meta_Results<-read.csv("metaOutputFDR_orderedByPval.csv", header=TRUE, stringsAsFactors = FALSE)

str(Meta_Results)
#'data.frame': 13814 obs. of  24 variables:

#Take a peek at your column names - for this code you will need the names for the Log2 Fold Change, p-value, and FDR columns for your meta-analysis results or variable of interest.

colnames(Meta_Results)

#############

#Read in the Volcano Plot Function:

VolcanoPlot<-function(DE_Results, VariableOfInterest, CoefficientCol, PvalueCol, FDRCol){

  #Grabbing the relevant columns of results and making a new dataframe with standardized column names:
  Coefficient<-DE_Results[[CoefficientCol]]
  Pvalue<-DE_Results[[PvalueCol]]
  FDR<-DE_Results[[FDRCol]]
  
  DE_Results_DF<-cbind.data.frame(Coefficient, Pvalue, FDR)  
  
  #Opening up a .tiff file to store the volcano plot
  tiff(paste("VolcanoPlot_", VariableOfInterest, ".tiff", sep=""), width = 5, height = 5, units = 'in', res = 300, compression = "lzw")
  
  par(mai=c(1.02, 1,0.9,0.40))

  #Making the volcano plot:

  with(DE_Results_DF, plot(Coefficient, -log10(Pvalue), pch=19, main=paste("Effect of", VariableOfInterest, sep=" "), cex.lab=1.8, cex.main=2, cex=0.6, xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
  
  #I've taken the parameters out for setting limits for the x and y axes:
  #xlim=c(-4.5,4.5), ylim=c(0,8),
  
  #Labeling different subsets of genes in the volcano plot with colors depending on their Log2FC, FDR, or both: 
  
  with(subset(DE_Results_DF, abs(Coefficient)>1), points(Coefficient, -log10(Pvalue), pch=19, col="red", cex=0.6))
  
  with(subset(DE_Results_DF, FDR<0.05), points(Coefficient, -log10(Pvalue), pch=19, col="green", cex=0.6))
  
  with(subset(DE_Results_DF, abs(Coefficient)>1 & FDR<0.05), points(Coefficient, -log10(Pvalue), pch=19, col="gold", cex=0.6))
  
  #I've taken the legend out of this function because it is hard to automate its placement and easy enough to write in the figure legend.
  #legend(-1.5, 7.6, legend=c("Log2FC > 1", "FDR<0.05", "both"), col=c("red", "green", "gold"), pch=19, cex=1)

  #closing the .tiff file connection so that the volcano plot is complete:
  dev.off()

}

###############

#Example Usage:

#Object containing differential expression or meta-analysis results:
DE_Results<-Meta_Results

#Categorical Variable of Interest:
VariableOfInterest<-"ADHD_Models"

#Name of column containing Log2 Fold Change for the Variable of Interest:
#Coefficient<-Meta_Results$Log2FC_estimate
CoefficientCol<-"Log2FC_estimate"

#Name of column containing p-value for the Variable of Interest:
PvalueCol<-"pval"

#Name of column containing FDR for the Variable of Interest:
FDRCol<-"FDR"

VolcanoPlot(DE_Results, VariableOfInterest, CoefficientCol, PvalueCol, FDRCol)
