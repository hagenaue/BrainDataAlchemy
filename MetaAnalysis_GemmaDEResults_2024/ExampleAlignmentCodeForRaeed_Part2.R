#We want to join this ortholog database to our mouse results (Log2FC and SV):

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(HOM_MouseVsRat, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")

#Structure of the new object:
str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:

Mouse_MetaAnalysis_SV_wOrthologs<-join(HOM_MouseVsRat, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")

#Structure of the new object:
str(Mouse_MetaAnalysis_SV_wOrthologs)
#'data.frame':	25288 obs. of  6 variables:

########################

#*If there are rat datasets*, we then want to join our mouse Log2FC and SV results to the rat results using the ortholog information:

MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")

#Structure of the new object:
str(MetaAnalysis_FoldChanges)
#'data.frame':	28101 obs. of  7 variables:

MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")

#Structure of the new object:
str(MetaAnalysis_SV)
#'data.frame':	28101 obs. of  7 variables:

#################

#Adding in the new GEO2R aligned rat data

#Rename the colnames so that the gene symbol column matches the name used in HOM_MouseVsRat:

colnames(NewRat_MetaAnalysis_FoldChanges)[1]<-"Rat_Symbol"
colnames(NewRat_MetaAnalysis_SV)[1]<-"Rat_Symbol"


MetaAnalysis_FoldChanges2<-join(MetaAnalysis_FoldChanges, NewRat_MetaAnalysis_FoldChanges, by="Rat_Symbol", type="full")

#Structure of the new object:
str(MetaAnalysis_FoldChanges2)

MetaAnalysis_SV2<-join(MetaAnalysis_SV, NewRat_MetaAnalysis_SV, by="Rat_Symbol", type="full")

#Structure of the new object:
str(MetaAnalysis_SV2)


#########################

#Once we confirm that code worked ok, we can rename the object to match the downstream code

MetaAnalysis_FoldChanges<-MetaAnalysis_FoldChanges2
MetaAnalysis_SV<-MetaAnalysis_SV2

#########################

#For simplicity's sake, I'm going to replace that Mouse-Rat Entrez annotation
#Because it is missing entries for any genes in the datasets that *don't* have orthologs
MetaAnalysis_FoldChanges$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_FoldChanges$Mouse_EntrezGene.ID, MetaAnalysis_FoldChanges$Rat_EntrezGene.ID, sep="_")

MetaAnalysis_SV$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_SV$Mouse_EntrezGene.ID, MetaAnalysis_SV$Rat_EntrezGene.ID, sep="_")

##############

#So currently we have more than 3 columns of annotation
#The original code is going to assume that there are only 3:
#RatEntrezID, MouseEntrezID, RatVsMouseEntrezID

#So we need to remove extra columns of annotation

colnames(MetaAnalysis_FoldChanges)

#Example: If columns 3:8 and 10:12 are all other types of gene annotation
#We could cut them out using code like this:
MetaAnalysis_FoldChanges<-MetaAnalysis_FoldChanges[,-c(3:8, 10:12)]
MetaAnalysis_SV<-MetaAnalysis_SV[,-c(3:8, 10:12)]
#But change the column numbers in this code to match what you see in the object

########################

###And then I think the rest of the meta-analysis code will just work.

#For the meta-analysis code the NumberOfComparisons and CutOffForNAs will be higher because you have more input now

