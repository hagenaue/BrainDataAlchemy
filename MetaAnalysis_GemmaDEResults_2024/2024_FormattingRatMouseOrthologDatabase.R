#Determining which mouse genes align to which rat genes
#This code was updated April 25, 2024

#Gene Ortholog Database:
#Downloaded from Jackson Labs April 25, 2024
#https://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt

#Another good source of ortholog information is HCOP:
#https://www.genenames.org/tools/hcop/

#I should note: when determining which genes align across species you will see two terms: Ortholog and Homolog
#"While homologous genes can be similar in sequence, similar sequences are not necessarily homologous. 
#Orthologous are homologous genes where a gene diverges after a speciation event, but the gene and its main function are conserved."
  
#Setting working directory:
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Teaching/BrainDataAlchemy/Summer_2024/Summer2024_Pipeline")

#Reading in the ortholog database:
HOM_AllOrganism<-read.delim("HOM_AllOrganism_20240425.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

colnames(HOM_AllOrganism)
# [1] "DB.Class.Key"                                       
# [2] "Common.Organism.Name"                               
# [3] "NCBI.Taxon.ID"                                      
# [4] "Symbol"                                             
# [5] "EntrezGene.ID"                                      
# [6] "Mouse.MGI.ID"                                       
# [7] "HGNC.ID"                                            
# [8] "OMIM.Gene.ID"                                       
# [9] "Genetic.Location"                                   
# [10] "Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7."
# [11] "Name"                                               
# [12] "Synonyms" 


table(HOM_AllOrganism$Common.Organism.Name)
# human mouse, laboratory               rat         zebrafish 
# 24586             21814             22466             29863 

HOM_Mouse<-HOM_AllOrganism[HOM_AllOrganism$Common.Organism.Name=="mouse, laboratory",]
str(HOM_Mouse)
#'data.frame':	21814 obs. of  12 variables:

HOM_Rat<-HOM_AllOrganism[HOM_AllOrganism$Common.Organism.Name=="rat",]
str(HOM_Rat)
#'data.frame':	22466 obs. of  12 variables:

sum(HOM_Mouse$DB.Class.Key%in%HOM_Rat$DB.Class.Key)
#[1] 19905
#So most of the mouse genes have orthologs in rats

#The column names are mostly the same in the two databases
#If I join them like that, R will get confused
colnames(HOM_Mouse)<-paste("Mouse", colnames(HOM_Mouse), sep="_")
colnames(HOM_Rat)<-paste("Rat", colnames(HOM_Rat), sep="_")

#But I still need the "shared" column name to be the same for the joining process
colnames(HOM_Mouse)[1]<-"DB.Class.Key"
colnames(HOM_Rat)[1]<-"DB.Class.Key"

library(plyr)
HOM_MouseVsRat<-join(HOM_Mouse, HOM_Rat, by="DB.Class.Key", type="full", match="all")
str(HOM_MouseVsRat)
#'data.frame':	24376 obs. of  23 variables:
#This data frame is bigger than the data frame for either rats or mice
#So some genes must either not have an ortholog or map to multiple genes in the other species

head(HOM_MouseVsRat)
colnames(HOM_MouseVsRat)

#SanityCheck
#Most (~75%) mouse and rat gene symbols should be the same
head(HOM_MouseVsRat[,c(4,15)])
# Mouse_Symbol Rat_Symbol
# 1        Wdr53      Wdr53
# 2          Tdg        Tdg
# 3       Trarg1     Trarg1
# 4        Pdgfb      Pdgfb
# 5       Gpr171     Gpr171
# 6      Glyatl3    Glyatl3

sum(HOM_MouseVsRat[,4]==HOM_MouseVsRat[,15], na.rm=TRUE)
#[1] 17288

#Looks good

write.csv(HOM_MouseVsRat, "HOM_MouseVsRat_20240425.csv")

#A quick check as to whether there is a 1-to-1 mapping between Entrez Gene ID and Gene Symbol
length(unique(HOM_MouseVsRat$Mouse_EntrezGene.ID))
#[1] 21813
length(unique(paste(HOM_MouseVsRat$Mouse_EntrezGene.ID, HOM_MouseVsRat$Mouse_Symbol, sep="_")))
#[1] 21813

length(unique(HOM_MouseVsRat$Rat_EntrezGene.ID))
#[1] 19184
length(unique(paste(HOM_MouseVsRat$Rat_EntrezGene.ID, HOM_MouseVsRat$Rat_Symbol, sep="_")))
#[1] 19184

#Are there repeats/multimaps?

max(table(HOM_MouseVsRat$Mouse_EntrezGene.ID))
#[1] 13

#Similar for symbol:
max(table(HOM_MouseVsRat[,4]))
#[1] 13
#Yes, at least one gene has 13 entries
tail(table(HOM_MouseVsRat[,4])[order(table(HOM_MouseVsRat[,4]))])
#Klk1b24 Klk1b26  Klk1b3  Klk1b5  Klk1b9  Klra17 
#11      11      11      11      11      13 

#how many genes are multimapped between species?
str(table(HOM_MouseVsRat$Mouse_EntrezGene.ID)[table(HOM_MouseVsRat$Mouse_EntrezGene.ID)>1])
# 'table' int [1:1299(1d)] 2 2 2 2 3 3 3 2 2 2 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:1299] "11465" "11837" "11951" "12268" ...

str(table(HOM_MouseVsRat$Mouse_EntrezGene.ID)[table(HOM_MouseVsRat$Mouse_EntrezGene.ID)>2])
# 'table' int [1:632(1d)] 3 3 3 3 4 3 5 11 10 11 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:632] "12313" "12314" "12315" "12491" ...

#632 mouse genes map to more than 2 rat genes

max(table(HOM_MouseVsRat$Rat_EntrezGene.ID))
#[1] 106

#same for gene symbol:
max(table(HOM_MouseVsRat[,15]))
#[1] 106
#Worse for rat
tail(table(HOM_MouseVsRat[,15])[order(table(HOM_MouseVsRat[,15]))])
# LOC685668    MGC116197 LOC100910078 LOC103694730 LOC108353098   Ssty1 
# 81           81          106          106          106          106 

#If a gene has 106 orthologs... I'm not so convinced that orthology is very meaningful.
#Let's get rid of them.

str(table(HOM_MouseVsRat$Rat_EntrezGene.ID)[table(HOM_MouseVsRat$Rat_EntrezGene.ID)>2])
# 'table' int [1:393(1d)] 3 3 11 3 3 4 4 4 4 15 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:393] "24242" "24244" "24256" "24282" ...
#393 rat genes map to more than 2 mouse gene

#Is it possible some of these are just NAs?
sum(is.na(HOM_MouseVsRat$Rat_EntrezGene.ID))
#[1] 1908

sum(is.na(HOM_MouseVsRat$Mouse_EntrezGene.ID))
#[1] 2

#Those aren't very helpful for identifying orthologs either...

#Let's get rid of the NAs first

HOM_MouseVsRat_noNA<-HOM_MouseVsRat[is.na(HOM_MouseVsRat$Rat_EntrezGene.ID)==FALSE & is.na(HOM_MouseVsRat$Mouse_EntrezGene.ID)==FALSE,]
str(HOM_MouseVsRat_noNA)
#'data.frame':	22466 obs. of  23 variables:

head(names(table(HOM_MouseVsRat$Mouse_EntrezGene.ID)[table(HOM_MouseVsRat$Mouse_EntrezGene.ID)>2]))
head(names(table(HOM_MouseVsRat$Rat_EntrezGene.ID)[table(HOM_MouseVsRat$Rat_EntrezGene.ID)>2]))

#These are the genes that I would like to keep:
MouseGenes_MappedLessThan2<-names(table(HOM_MouseVsRat$Mouse_EntrezGene.ID)[table(HOM_MouseVsRat$Mouse_EntrezGene.ID)<2])
str(MouseGenes_MappedLessThan2)
#chr [1:20530] 

RatGenes_MappedLessThan2<-names(table(HOM_MouseVsRat$Rat_EntrezGene.ID)[table(HOM_MouseVsRat$Rat_EntrezGene.ID)<2])
str(RatGenes_MappedLessThan2) 
#chr [1:18163] 
 
#The fact that EntrezID is encoded as a character here but as an integer in the HOM_MouseVsRat is going to cause issues
#This will come up in other situations as well

#I may just force EntrezID into a character format:
HOM_MouseVsRat_noNA$Mouse_EntrezGene.ID<-as.character(HOM_MouseVsRat_noNA$Mouse_EntrezGene.ID)
HOM_MouseVsRat_noNA$Rat_EntrezGene.ID<-as.character(HOM_MouseVsRat_noNA$Rat_EntrezGene.ID)
str(HOM_MouseVsRat_noNA)
#Looks good

HOM_MouseVsRat_noNA_LessThan2Orthologs<-HOM_MouseVsRat_noNA[(HOM_MouseVsRat_noNA$Mouse_EntrezGene.ID%in%MouseGenes_MappedLessThan2) & (HOM_MouseVsRat_noNA$Rat_EntrezGene.ID%in%RatGenes_MappedLessThan2),]
str(HOM_MouseVsRat_noNA_LessThan2Orthologs)
#'data.frame':	17241 obs. of  23 variables:

#Out of curiosity, how much have I gained by switching to this method vs. just matching gene symbol?
sum(HOM_MouseVsRat_noNA_LessThan2Orthologs$Mouse_Symbol==HOM_MouseVsRat_noNA_LessThan2Orthologs$Rat_Symbol)
#[1] 16479
16479/17241
#[1] 0.955803
#So using the simpler method of just matching gene symbol got us 95% to the same point. Sigh.

HOM_MouseVsRat_noNA_LessThan2Orthologs$MouseRat_EntrezGene.ID<-paste(HOM_MouseVsRat_noNA_LessThan2Orthologs$Mouse_EntrezGene.ID, HOM_MouseVsRat_noNA_LessThan2Orthologs$Rat_EntrezGene.ID, sep="_")

MouseVsRat_NCBI_Entrez<-cbind.data.frame(Rat_EntrezGene.ID=HOM_MouseVsRat_noNA_LessThan2Orthologs$Rat_EntrezGene.ID, Mouse_EntrezGene.ID=HOM_MouseVsRat_noNA_LessThan2Orthologs$Mouse_EntrezGene.ID, MouseVsRat_EntrezGene.ID=HOM_MouseVsRat_noNA_LessThan2Orthologs$MouseRat_EntrezGene.ID)

str(MouseVsRat_NCBI_Entrez)
# 'data.frame':	17246 obs. of  3 variables:
# $ Rat_EntrezGene.ID       : chr  "498097" "114521" "360576" "24628" ...
# $ Mouse_EntrezGene.ID     : chr  "68980" "21665" "237858" "18591" ...
# $ MouseVsRat_EntrezGene.ID: chr  "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...

tail(MouseVsRat_NCBI_Entrez)

sum(is.na(MouseVsRat_NCBI_Entrez[,3]))

write.csv(MouseVsRat_NCBI_Entrez, "MouseVsRat_NCBI_Entrez_JacksonLab_20240425.csv")


