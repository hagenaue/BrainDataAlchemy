#Determining which mouse genes align to which rat genes
#This code was updated May 11, 2026

#Gene Ortholog Database:
#Downloaded from Jackson Labs May 11, 2026
#https://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt
#https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt 

#Another good source of ortholog information is HCOP:
#https://www.genenames.org/tools/hcop/

#I should note: when determining which genes align across species you will see two terms: Ortholog and Homolog
#"While homologous genes can be similar in sequence, similar sequences are not necessarily homologous. 
#Orthologous are homologous genes where a gene diverges after a speciation event, but the gene and its main function are conserved."

#Setting working directory:
setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/CodeSharing")

#Reading in the ortholog database:
HOM_AllOrganism<-read.delim("HOM_AllOrganism_20260511.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

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
# human  mouse, laboratory            rat         zebrafish 
# 24592             21930             23651             26036 

HOM_Mouse<-HOM_AllOrganism[HOM_AllOrganism$Common.Organism.Name=="mouse, laboratory",]
str(HOM_Mouse)
#'data.frame':	21930 obs. of  12 variables:

HOM_Rat<-HOM_AllOrganism[HOM_AllOrganism$Common.Organism.Name=="rat",]
str(HOM_Rat)
#'data.frame':	23651 obs. of  12 variables:

sum(HOM_Mouse$DB.Class.Key%in%HOM_Rat$DB.Class.Key)
#[1] 20123
#So most of the mouse genes have orthologs in rats

#It would make our life easier for other comparisons if we also had easy access to the human orthologs:
HOM_Human<-HOM_AllOrganism[HOM_AllOrganism$Common.Organism.Name=="human",]
str(HOM_Human)
#'data.frame':	24592 obs. of  12 variables:

sum(HOM_Mouse$DB.Class.Key%in%HOM_Human$DB.Class.Key)
#[1] 20184

#The column names are mostly the same in the two databases
#If I join them like that, R will get confused
colnames(HOM_Mouse)<-paste("Mouse", colnames(HOM_Mouse), sep="_")
colnames(HOM_Rat)<-paste("Rat", colnames(HOM_Rat), sep="_")
colnames(HOM_Human)<-paste("Human", colnames(HOM_Human), sep="_")

#But I still need the "shared" column name to be the same for the joining process
colnames(HOM_Mouse)[1]<-"DB.Class.Key"
colnames(HOM_Rat)[1]<-"DB.Class.Key"
colnames(HOM_Human)[1]<-"DB.Class.Key"

head(HOM_Mouse[,1])

colnames(HOM_Mouse)

#"All genes sharing the same DB.Class.Key number are considered part of the same orthology cluster"
head(HOM_Mouse$DB.Class.Key)
#[1] 51424856 51424857 51424858 51424859 51424860 51424861
sum(is.na(HOM_Mouse$DB.Class.Key))
#[1] 0
table(table(HOM_Mouse$DB.Class.Key))
# 1     2 
# 21924     3 

head(HOM_Mouse$Mouse_Mouse.MGI.ID)
#[1] "MGI:1340024" "MGI:98660"   "MGI:98360"   "MGI:2443895" "MGI:98001"   "MGI:894689" 
sum(is.na(HOM_Mouse$Mouse_Mouse.MGI.ID))
#[1] 0
table(table(HOM_Mouse$Mouse_Mouse.MGI.ID))
# 1     2 
# 21924     3 

##################

#Adding in Ensembl annotation (it seems to be missing)

MRK_ENSEMBL<-read.delim("MRK_ENSEMBL_20260511.txt", header=FALSE, stringsAsFactors = FALSE)
colnames(MRK_ENSEMBL)
str(MRK_ENSEMBL)

# 'data.frame':	77645 obs. of  13 variables:
# $ V1 : chr  "MGI:1915733" "MGI:1919275" "MGI:1914753" "MGI:1916606" ...
# $ V2 : chr  "1110002O04Rik" "1600012P17Rik" "1700001G17Rik" "1700003I22Rik" ...
# $ V3 : chr  "RIKEN cDNA 1110002O04 gene" "RIKEN cDNA 1600012P17 gene" "RIKEN cDNA 1700001G17 gene" "RIKEN cDNA 1700003I22 gene" ...
# $ V4 : num  -1 68.5 12.8 -1 -1 ...
# $ V5 : chr  "1" "1" "1" "1" ...
# $ V6 : chr  "ENSMUSG00000102531" "ENSMUSG00000047661" "ENSMUSG00000103746" "ENSMUSG00000100372" ...
# $ V7 : chr  "ENSMUST00000268248 ENSMUST00000268252 ENSMUST00000194261 ENSMUST00000268253 ENSMUST00000268247 ENSMUST000002682"| __truncated__ "ENSMUST00000062159 ENSMUST00000162474" "" "ENSMUST00000332611 ENSMUST00000190280 ENSMUST00000186048" ...
# $ V8 : chr  "" "" "" "" ...
# $ V9 : chr  "lncRNA gene" "lncRNA gene" "lncRNA gene" "lncRNA gene" ...
# $ V10: int  35909029 158795271 33708905 56058137 137253172 186857335 120363513 80423649 53197736 138804010 ...
# $ V11: int  35920244 158808032 33709793 56059362 137253580 186860049 120393663 80453377 53226795 138804688 ...
# $ V12: chr  "+" "-" "+" "+" ...
# $ V13: chr  "ncRNA|lncRNA" "lncRNA|ncRNA" "TEC|ncRNA" "ncRNA|lncRNA" ...

#This only provides us with the ENSEMBL gene annotation for mice :(

table(table(MRK_ENSEMBL$V1))
# 1     2     3     5 
# 77312   161     2     1 

table(table(MRK_ENSEMBL$V2))
# 1     2     3     5 
# 77312   161     2     1 

table(table(MRK_ENSEMBL$V6))
#    1 
# 77645 

#Are there really 77,645 Ensembl mouse genes??
#Yep - double-checked and there are.
#And apparently they all have symbols???

MRK_ENSEMBL$V2[c(1:500)]
#Looks like the ambiguous ones have Riken symbols.


#Also, this is a lot of info - let's trim it down:

MRK_ENSEMBL2<-MRK_ENSEMBL[,c(1,2,3,6)]

colnames(MRK_ENSEMBL2)<-c("Mouse_Mouse.MGI.ID", "Mouse_Symbol", "Mouse_Name", "Mouse_ENSEMBLGene.ID")

dim(MRK_ENSEMBL2)
#[1] 77645     4

sum(is.na(MRK_ENSEMBL2$Mouse_Mouse.MGI.ID))
#[1] 0

table(table(MRK_ENSEMBL2$Mouse_Mouse.MGI.ID))
# 1        2     3     5 
# 77312   161     2     1 

sum(is.na(MRK_ENSEMBL2$Mouse_ENSEMBLGene.ID))
#[1] 0

sum(is.na(MRK_ENSEMBL2$Mouse_Symbol))
#[1] 0

sum(is.na(HOM_Mouse))

HOM_Mouse2<-join(HOM_Mouse, MRK_ENSEMBL2, by="Mouse_Mouse.MGI.ID", type="full")
str(HOM_Mouse2)
#'data.frame':	77713 obs. of  13 variables:

sum(is.na(HOM_Mouse2$Mouse_Mouse.MGI.ID))
#[1] 0

sum(is.na(HOM_Mouse2$Mouse_Symbol))
#[1] 0

table(table(HOM_Mouse2$Mouse_Symbol))
#    1     2     3     5 
# 77374   164     2     1 

sum(is.na(HOM_Mouse2$Mouse_EntrezGene.ID))
#[1] 55749

sum(is.na(HOM_Mouse2$Mouse_ENSEMBLGene.ID))
#[1] 66

table(table(HOM_Mouse2$Mouse_ENSEMBLGene.ID))
# 1     2 
# 77643     2

sum(is.na(HOM_Mouse2$Mouse_EntrezGene.ID))
#[1] 55749

##################

#Grabbing only the mouse genes that both Ensembl and Entrez agree exist:

HOM_Mouse2_NoNA<-HOM_Mouse2[((is.na(HOM_Mouse2$Mouse_EntrezGene.ID)==FALSE) & (is.na(HOM_Mouse2$Mouse_ENSEMBLGene.ID)==FALSE)),]

dim(HOM_Mouse2_NoNA)
#[1] 21898    14

HOM_Mouse2_NoNA$Ensembl_Entrez<-paste(HOM_Mouse2_NoNA$Mouse_ENSEMBLGene.ID, HOM_Mouse2_NoNA$Mouse_EntrezGene.ID, sep="_")
  
table(table(HOM_Mouse2_NoNA$Ensembl_Entrez))
# 1         2 
# 21896     1 

table(table(HOM_Mouse2_NoNA$Mouse_EntrezGene.ID))
#  1        2     3 
# 21830    31     2

#Removing any genes that have Entrez mapping to multiple Ensembl or vice-versa:

HOM_Mouse2_NoNA_Clean<-HOM_Mouse2_NoNA[((duplicated(HOM_Mouse2_NoNA$Mouse_EntrezGene.ID)==FALSE) & (duplicated(HOM_Mouse2_NoNA$Mouse_ENSEMBLGene.ID)==FALSE)),]

table(table(HOM_Mouse2_NoNA_Clean$Mouse_EntrezGene.ID))
#    1 
#21862 

table(table(HOM_Mouse2_NoNA_Clean$Mouse_ENSEMBLGene.ID))
#     1 
# 21862 

colnames(HOM_Mouse2_NoNA_Clean)

write.csv(HOM_Mouse2_NoNA_Clean, "Mouse_GeneAnnotation_Ensembl_Entrez_Agreed.csv")

##################

#Adding Rat Ensembl annotation to the data:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Rn.eg.db")


library(org.Rn.eg.db)


#Outputting gene symbols for all ENSEMBL genes:

ls("package:org.Rn.eg.db")

#grabbing all of the ensembl ids in the annotation package:
uniKeys<-keys(org.Rn.eg.db, keytype="ENSEMBL")
str(uniKeys)
#chr [1:20998] "ENSRNOG00000028896" "ENSRNOG00000032908" "ENSRNOG00000009845" "ENSRNOG00000016924"...
#This is a vector with all ensembl genes in the org.Rn.eg.db

uniKeys[1]
# [1] "ENSRNOG00000028896"

sum(uniKeys == "ENSRNOG00000028896")


#We would like the symbol for all of those ensembl ids
cols<-c("SYMBOL", "GENENAME", "ENTREZID", "ENSEMBL")  

#grabbing symbol for every ensembl id in the package:
EnsemblVsGeneSymbol<-AnnotationDbi::select(org.Rn.eg.db, keys=uniKeys, columns=cols, keytype="ENSEMBL")  
dim(EnsemblVsGeneSymbol)
#[1] 40304     4

head(EnsemblVsGeneSymbol)
#              ENSEMBL SYMBOL                            GENENAME ENTREZID
# 1 ENSRNOG00000028896    A2m               alpha-2-macroglobulin    24153
# 2 ENSRNOG00000032908 Acaa1a       acetyl-CoA acyltransferase 1A    24157
# 3 ENSRNOG00000009845  Acadm acyl-CoA dehydrogenase medium chain    24158
# 4 ENSRNOG00000016924   Acly                   ATP citrate lyase    24159
# 5 ENSRNOG00000005260   Acp1                  acid phosphatase 1    24161
# 6 ENSRNOG00000013594   Acp2       acid phosphatase 2, lysosomal    24162

str(EnsemblVsGeneSymbol)
# 'data.frame':	40304 obs. of  4 variables:
#   $ ENSEMBL : chr  "ENSRNOG00000028896" "ENSRNOG00000032908" "ENSRNOG00000009845" "ENSRNOG00000016924" ...
# $ SYMBOL  : chr  "A2m" "Acaa1a" "Acadm" "Acly" ...
# $ GENENAME: chr  "alpha-2-macroglobulin" "acetyl-CoA acyltransferase 1A" "acyl-CoA dehydrogenase medium chain" "ATP citrate lyase" ...
# $ ENTREZID: chr  "24153" "24157" "24158" "24159" ...

sum(is.na(EnsemblVsGeneSymbol$ENTREZID))
#[1] 0

table(table(EnsemblVsGeneSymbol$ENTREZID))
#     1     2     3     4     5     6     8    27    49    71    91 
# 25666    79    20     9    17     6     2    27     4    71    91

sum(is.na(EnsemblVsGeneSymbol$ENSEMBL))
#[1] 0

table(table(EnsemblVsGeneSymbol$ENSEMBL))
#   1     2     3     4     5     6     7     8    10    11    13    15    27    71    91 
# 25102   259    41    24    60     7     1     2     1     1     1     1    27    71    91 
table(table(EnsemblVsGeneSymbol$SYMBOL))
#     1     2     3     4     5     6     8    27    49    71    91 
# 25666    79    20     9    17     6     2    27     4    71    91 


#It seems...fishy... that there is a similar distribution of some of the high numbers - e.g. all 3 annotations have 71 genes that are repeated 71 times and 91 genes that are repeated 91 times.

table(EnsemblVsGeneSymbol$SYMBOL)[table(EnsemblVsGeneSymbol$SYMBOL)==91]
#LOC120100240 LOC120100241 LOC120100242 LOC120100243 LOC120100244 LOC120100246 LOC120100247 #...etc

EnsemblVsGeneSymbol[EnsemblVsGeneSymbol$SYMBOL=="LOC120100240",]
# 23749 ENSRNOG00000088697 LOC120100240 small nucleolar RNA SNORD115 120100240
# 23840 ENSRNOG00000078089 LOC120100240 small nucleolar RNA SNORD115 120100240
# 23931 ENSRNOG00000088825 LOC120100240 small nucleolar RNA SNORD115 120100240

table(EnsemblVsGeneSymbol$ENSEMBL)[table(EnsemblVsGeneSymbol$ENSEMBL)==91]
#ENSRNOG00000071337 ENSRNOG00000071680 ENSRNOG00000071840 ENSRNOG00000071997...

EnsemblVsGeneSymbol[EnsemblVsGeneSymbol$ENSEMBL=="ENSRNOG00000071337",]
# 24841 ENSRNOG00000071337 LOC120100240 small nucleolar RNA SNORD115 120100240
# 24842 ENSRNOG00000071337 LOC120100241 small nucleolar RNA SNORD115 120100241
# 24843 ENSRNOG00000071337 LOC120100242 small nucleolar RNA SNORD115 120100242

table(EnsemblVsGeneSymbol$ENTREZID)[table(EnsemblVsGeneSymbol$ENTREZID)==91]
#120100240 120100241 120100242 120100243 120100244 120100246 120100247 120100248 120100250...

head(duplicated(EnsemblVsGeneSymbol))

sum(duplicated(EnsemblVsGeneSymbol))
#[1] 0

#Grrrr...
#It looks like there are a lot of Ensembl and Entrez annotations that are in disagreement.

##################

#Joining the homology database with the rat Ensembl/Entrez mapping

dim(HOM_Rat)
#[1] 23651    12

colnames(HOM_Rat)

colnames(EnsemblVsGeneSymbol)<-c("Rat_Ensembl", "Rat_DBSymbol", "Rat_DBName", "Rat_EntrezGene.ID")

HOM_Rat2<-join(HOM_Rat, EnsemblVsGeneSymbol, by="Rat_EntrezGene.ID", type="left", match="all")

dim(HOM_Rat2)
#[1] 23880    15

#Grabbing only the mouse genes that both Ensembl and Entrez agree exist:

HOM_Rat_NoNA<-HOM_Rat2[((is.na(HOM_Rat2$Rat_EntrezGene.ID)==FALSE) & (is.na(HOM_Rat2$Rat_Ensembl)==FALSE)),]

dim(HOM_Rat_NoNA)
#[1] 21934    15

HOM_Rat_NoNA$Ensembl_Entrez<-paste(HOM_Rat_NoNA$Rat_Ensembl, HOM_Rat_NoNA$Rat_EntrezGene.ID, sep="_")

table(table(HOM_Rat_NoNA$Rat_EntrezGene.ID))
# 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    21 
# 18062   539   156    66    59    25    10    14     9     6     4     8     6     2     7     1     6     1 
# 23    25    27    28    43    60    65    76    77    97 
# 1     6     1     2     1     1     3     1     1     1 

#Removing any genes that have Entrez mapping to multiple Ensembl or vice-versa:

HOM_Rat_NoNA_Clean<-HOM_Rat_NoNA[((duplicated(HOM_Rat_NoNA$Rat_Ensembl)==FALSE) & (duplicated(HOM_Rat_NoNA$Rat_EntrezGene.ID)==FALSE)),]

table(table(HOM_Rat_NoNA_Clean$Rat_EntrezGene.ID))
#      1 
# 18860 

table(table(HOM_Rat_NoNA_Clean$Rat_Ensembl))
#    1 
# 18860   

colnames(HOM_Rat_NoNA_Clean)

write.csv(HOM_Rat_NoNA_Clean, "Rat_GeneAnnotation_Ensembl_Entrez_Agreed.csv")

###################

#Making a Rat vs Mouse ortholog database (version 1 - No Ensembl annotation)

library(plyr)

HOM_MouseVsRat<-join(HOM_Mouse, HOM_Rat, by="DB.Class.Key", type="full", match="all")

str(HOM_MouseVsRat)
#'data.frame':	25461 obs. of  23 variables:
#This data frame is bigger than the data frame for either rats or mice
#So some genes must either not have an ortholog or map to multiple genes in the other species

head(HOM_MouseVsRat)
colnames(HOM_MouseVsRat)

#SanityCheck
#Most (~75%) mouse and rat gene symbols should be the same
head(HOM_MouseVsRat[,c(4,15)])
# Mouse_Symbol Rat_Symbol
# Mouse_Symbol Rat_Symbol
# 1      Aldh1l1    Aldh1l1
# 2          Sry        Sry
# 3          Sry      Sry3A
# 4        Sox12      Sox12
# 5          Fry        Fry
# 6        Rpe65      Rpe65

sum(HOM_MouseVsRat[,4]==HOM_MouseVsRat[,15], na.rm=TRUE)
#[1] 17294

#Looks good

write.csv(HOM_MouseVsRat, "HOM_MouseVsRat_20260511.csv")

##################

#Making a Rat vs Mouse ortholog database (version 2 - with Ensembl annotation)

library(plyr)

HOM_MouseVsRat_NoNA_Clean<-join(HOM_Mouse2_NoNA_Clean, HOM_Rat_NoNA_Clean, by="DB.Class.Key", type="full", match="all")

str(HOM_MouseVsRat_NoNA_Clean)
#'data.frame':	22695 obs. of  28 variables:
#This data frame is bigger than the data frame for either rats or mice
#So some genes must either not have an ortholog or map to multiple genes in the other species

head(HOM_MouseVsRat_NoNA_Clean)
colnames(HOM_MouseVsRat_NoNA_Clean)

#SanityCheck
#Most (~75%) mouse and rat gene symbols should be the same
head(HOM_MouseVsRat_NoNA_Clean[,c(4,17)])
# Mouse_Symbol Rat_Symbol
# 1      Aldh1l1    Aldh1l1
# 2          Sry        Sry
# 3          Sry      Sry3A
# 4        Sox12      Sox12
# 5          Fry        Fry
# 6        Rpe65      Rpe65

sum(HOM_MouseVsRat_NoNA_Clean[,4]==HOM_MouseVsRat_NoNA_Clean[,17], na.rm=TRUE)
#[1] 16734

16734/22695
#[1] 0.737343
#Looks good

sum(is.na(HOM_MouseVsRat_NoNA_Clean$Mouse_EntrezGene.ID))
#[1] 22

sum(is.na(HOM_MouseVsRat_NoNA_Clean$Mouse_ENSEMBLGene.ID))
#[1] 22

sum(is.na(HOM_MouseVsRat_NoNA_Clean$Rat_Ensembl))
#[1] 3835

sum(is.na(HOM_MouseVsRat_NoNA_Clean$Rat_EntrezGene.ID))
#[1] 3835

3835+22
#[1] 3857

16734+3857
#[1] 20591

#So there are some multimapped genes between species.

write.csv(HOM_MouseVsRat_NoNA_Clean, "HOM_MouseVsRat_EntrezEnsemblAgree_20260511.csv")

########################

#There are some multimapped genes between species.

table(table(HOM_MouseVsRat_NoNA_Clean$Rat_EntrezGene.ID))
# 1 
# 18860 

table(table(HOM_MouseVsRat_NoNA_Clean$Mouse_EntrezGene.ID))
#   1     2     3     4     5     6     7     9    10    11    14 
# 21299   419    94    28    13     2     3     1     1     1     1 

HOM_MouseVsRat_NoNA_Clean_NoMultiMapped<-HOM_MouseVsRat_NoNA_Clean[((duplicated(HOM_MouseVsRat_NoNA_Clean$Mouse_EntrezGene.ID)==FALSE) | (duplicated(HOM_MouseVsRat_NoNA_Clean$Rat_EntrezGene.ID)==FALSE)),]

dim(HOM_MouseVsRat_NoNA_Clean_NoMultiMapped)
#[1] 22695    28

sum(HOM_MouseVsRat_NoNA_Clean_NoMultiMapped$Mouse_Symbol==HOM_MouseVsRat_NoNA_Clean_NoMultiMapped$Rat_Symbol, na.rm=TRUE)
#[1] 16734

write.csv(HOM_MouseVsRat_NoNA_Clean_NoMultiMapped, "HOM_MouseVsRat_EntrezEnsemblAgree_NoMultimapped_20260511.csv")

##################
##################
##################

#this code is leftover from May 11 2026 (before I cleaned up some of the code above):

HOM_MouseVsRatVsHuman<-join_all(list(HOM_Mouse, HOM_Rat, HOM_Human), by="DB.Class.Key", type="full", match="all")
str(HOM_MouseVsRatVsHuman)
#'data.frame':	33970 obs. of  34 variables:
#'Wow, a lot more non-matches or multimatches added in for humans.

colnames(HOM_MouseVsRatVsHuman)
head(HOM_MouseVsRatVsHuman[,c(4,15,26)])
# Mouse_Symbol Rat_Symbol Human_Symbol
# 1      Aldh1l1    Aldh1l1      ALDH1L1
# 2          Sry        Sry          SRY
# 3          Sry      Sry3A          SRY
# 4        Sox12      Sox12        SOX12
# 5          Fry        Fry          FRY
# 6        Rpe65      Rpe65        RPE65

sum(HOM_MouseVsRatVsHuman[,4]==HOM_MouseVsRatVsHuman[,15], na.rm=TRUE)
#[1] 18245
#So something like 1000 genes that have symbols in both mice and rats are multimatches with human genes.

write.csv(HOM_MouseVsRatVsHuman, "HOM_MouseVsRatVsHuman_20260511.csv")

##########

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

save.image("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/CodeSharing/Workspace_NewOrthologDB_20260521.RData")