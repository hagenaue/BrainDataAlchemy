#This code contains a function for making nice forest plots illustrating the effect sizes (Log2FC) for each of the statistical contrasts and datasets included in our meta-analysis for a particular gene (e.g., one of our top findings)
#Megan Hagenauer
#July 25 2024

##########################

library(metafor)

##########################

#Function:

MakeForestPlots<-function(metaOutputFDR_annotated, EntrezIDAsCharacter, species){
  
  #I originally wrote this function using only mouse Entrez IDs as input 
  #but then I realized that the function didn't work for genes that were only found in rats
  #so now the function allows either rat or mouse Entrez ids as input, and includes a conditional (if/else) statement
  
  if(species=="Mouse"){
    
    #This grabs the mouse gene symbol for the EntrezID from our meta-analysis annotation:
    MouseGeneSymbol<-metaOutputFDR_annotated$Mouse_Symbol[which(metaOutputFDR_annotated$Mouse_EntrezGene.ID==EntrezIDAsCharacter)]
    #This grabs the rat gene symbol for the EntrezID from our meta-analysis annotation:
    RatGeneSymbol<-metaOutputFDR_annotated$Rat_Symbol[which(metaOutputFDR_annotated$Mouse_EntrezGene.ID==EntrezIDAsCharacter)]
   
    #This grabs the Log2FC values for the EntrezID
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter),-c(1:3)])
    
    #This grabs the sampling variance (SV) values for the EntrezID
    var<-as.numeric(MetaAnalysis_SV_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter),-c(1:3)]) 
    
  }else if(species=="Rat"){
    
    #This set of code does all of the same processes as above, but interpreting the EntrezID as a Rat Entrez ID:
    
    RatGeneSymbol<-metaOutputFDR_annotated$Rat_Symbol[which(metaOutputFDR_annotated$Rat_EntrezGene.ID==EntrezIDAsCharacter)]
    
    MouseGeneSymbol<-metaOutputFDR_annotated$Mouse_Symbol[which(metaOutputFDR_annotated$Rat_EntrezGene.ID==EntrezIDAsCharacter)]
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Rat_EntrezGene.ID==EntrezIDAsCharacter),-c(1:3)])
    
    var<-as.numeric(MetaAnalysis_SV_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Rat_EntrezGene.ID==EntrezIDAsCharacter),-c(1:3)]) 
    
  }else{
    
    print("Please use either 'Mouse' or 'Rat' to indicate whether you are using annotation for mouse or rat genes")
    
  }
    
#This code makes the Forest Plot
  
  #First it opens up a .pdf file to output the plot into:
  #It automatically names that file with the mouse and rat ene symbols
  pdf(paste("ForestPlot_Mouse", MouseGeneSymbol,"Rat", RatGeneSymbol,".pdf", sep="_"), height=5, width=8)
  
  #This code makes the forest plot:
  #Note that the x-axis limits are currently set to -3 to 3
  #This may be too big or too small for visualizing the results for some genes.
  forest.rma(rma(effect, var), slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)],  xlim=c(-3, 3))
  
  #This code labels the forest plot with the mouse and rat gene symbols:
  mtext(paste("Mouse", MouseGeneSymbol, "Rat", RatGeneSymbol, sep="_"), line=-1.5, cex=2)
  
  #This closes the connection to the .pdf file, finishing the plot
  dev.off()
  
}

#Note: the way that this function is currently written, I suspect it might potentially throw up an error message if there is more than one set of log2FC associated with an Entrez ID. 
#This would happen if an Entrez ID in one species (e.g., mouse) mapped to more than one Entrez ID in the other species (e.g., rat) at the point that the results from the two species was joined. 
#It should be solvable by just using the EntrezID for the other species (e.g., rat) as the input to the function


#######################

#Example Usage:

#Note - this function currently uses Mouse Entrez ID (NCBI ID) as it's input
#It needs to be formatted as a character (not an integer) to work
#I will need to make a version of this later that accepts rat Entrez ID

#MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="13170", species="Mouse") #Dbp

#MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="74772", species="Mouse") #Atp13a2

#MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="282580", species="Rat") #Stab2




