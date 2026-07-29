#This code contains a function for making nice forest plots illustrating the effect sizes (Log2FC) for each of the statistical contrasts and datasets included in our meta-analysis for a particular gene (e.g., one of our top findings)
#Megan Hagenauer
#July 25 2024
#Updated July 29, 2026

##########################

if (!require("metafor", quietly = TRUE)){
  install.packages("metafor")
}

library(metafor)

##########################

#Function:

MakeForestPlots<-function(metaOutputFDR_annotated, Combined_MouseRat_GeneSymbol, Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches){
   
    #This grabs the Log2FC values for the Gene
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$MouseRat_GeneSymbol==Combined_MouseRat_GeneSymbol),Columns_DE])
    
    #This grabs the sampling variance (SV) values for the Gene
    var<-as.numeric(MetaAnalysis_SV_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$MouseRat_GeneSymbol==Combined_MouseRat_GeneSymbol),Columns_DE]) 
    
#This code makes the Forest Plot
  
  #First it opens up a .pdf file to output the plot into:
  #It automatically names that file with the mouse and rat ene symbols
  pdf(paste("ForestPlot", Combined_MouseRat_GeneSymbol,".pdf", sep="_"), height=HeightInInches, width=8)
  
  #This code makes the forest plot:
  #Note that the x-axis limits are currently set to -3 to 3
  #This may be too big or too small for visualizing the results for some genes.
  forest.rma(rma(effect, var), slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[Columns_DE], xlim=Xaxis_LowerAndUpperBound)
  
  #This code labels the forest plot with the mouse and rat gene symbols:
  mtext(paste(Combined_MouseRat_GeneSymbol, sep=""), line=-1.5, cex=2)
  
  #This closes the connection to the .pdf file, finishing the plot
  dev.off()
  
}



#######################

#Example Usage:

#Xaxis_LowerAndUpperBound<-c(-6,6)
#HeightInInches<-5

#MakeForestPlots(metaOutputFDR_annotated, "Frem3_Frem3", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

#MakeForestPlots(metaOutputFDR_annotated, "Scin_Scin", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)



