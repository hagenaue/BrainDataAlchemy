#This is R code for a practice search in GEO
#2026-06-22
#Megan Hagenauer

########################

if (!require(BiocManager, quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

#Come back to protocols.io and add quotes to "BiocManager" and "GEOquery"
#protocols.io: are quotation marks in code a problem?

########################

searchGEO

#Example search code:
MyQueryTerms<-'((environ*[All Fields] AND enrich*[All Fields]) OR (enrich*[All Fields] AND housing[All Fields]) OR (enrich*[All Fields] AND housed[All Fields]) OR (social*[All Fields] AND housing[All Fields]) OR (social*[All Fields] AND housed[All Fields]) OR "run"[All Fields] OR running[All Fields] OR exercis*[All Fields] OR wheel*[All Fields] OR toy*[All Fields] OR welfare[All Fields] OR (social[All Fields] AND enrich*[All Fields]) OR (sensor*[All Fields] AND enrich*[All Fields]) OR (motor[All Fields] AND enrich*[All Fields]) OR (cognitiv*[All Fields] AND enrich*[All Fields]) OR (behav*[All Fields] AND enrich*[All Fields]) OR (experienc*[All Fields] AND novel*[All Fields]) OR (environmen*[All Fields] AND novel*[All Fields]) OR (stimulat*[All Fields] AND novel*[All Fields]) OR (stimulat*[All Fields] AND environmen*[All Fields]) OR (stimulat*[All Fields] AND social*[All Fields]) OR (stimulat*[All Fields] AND cognitiv*[All Fields]) OR (stimulat*[All Fields] AND motor*[All Fields]) OR lifestyle[All Fields]) AND (hippocamp*[All Fields] OR "dentate gyrus"[All Fields] OR CA1[All Fields] OR CA2[All Fields] OR CA3[All Fields] OR "cornu ammonis"[All Fields]) AND ("Mus musculus"[ORGN] OR "Rattus norvegicus"[ORGN]) AND ("Expression profiling by high throughput sequencing"[DataSet Type] OR "Expression profiling by array"[DataSet Type]) AND "gse"[Filter]' 

#Another example search code:
MyQueryTerms<-'(((chronic*[All Fields] OR subchronic*[All Fields] OR “sub-chronic”[All Fields] OR repeated[All Fields] OR repetitive[All Fields] OR prolonged[All Fields] OR sustained[All Fields] OR continuous*[All Fields] OR extended[All Fields] OR distress*[All Fields]) AND (stress*[All Fields] OR defeat*[All Fields] OR restrain*[All Fields] OR immobil*[All Fields] OR advers*[All Fields] OR subordination[All Fields] OR isolation[All Fields])) OR (allostasis[All Fields]) OR ("allostatic load"[All Fields]) OR ("learned helplessness" [All Fields]) OR (PTSD[All Fields]) OR (sCVS[All Fields]) OR(CVS[All Fields]) OR (CSDS[All Fields]) OR (SDS[All Fields]) OR (“CUS”[All Fields]) OR (CUMS[All Fields]) OR (CMS[All Fields])) AND (hippocamp*[All Fields] OR "dentate gyrus"[All Fields] OR CA1[All Fields] OR CA2[All Fields] OR CA3[All Fields] OR "cornu ammonis"[All Fields]) AND ("Mus musculus"[ORGN] OR "Rattus norvegicus"[ORGN]) AND ("Expression profiling by high throughput sequencing"[DataSet Type] OR "Expression profiling by array"[DataSet Type]) AND "gse"[Filter]' 


QueryResults <- searchGEO(MyQueryTerms)

str(QueryResults)

# This will show you an overview of the identified GEO Records:

str(QueryResults)


#Adding columns to hold the additional metadata to our Query Results object:

QueryResults$Citation<-character(length=nrow(QueryResults))
QueryResults$PMID<-character(length=nrow(QueryResults))
QueryResults$Contributor<-character(length=nrow(QueryResults))
QueryResults$Date<-character(length=nrow(QueryResults))
QueryResults$Abstract<-character(length=nrow(QueryResults))

#Looping over each of the identified GEO records and extracting the desired metadata:


for(i in c(1:nrow(QueryResults))){
  
  gse_raw <- getGEO(QueryResults$`Series Accession`[i], GSEMatrix=FALSE)
  
  QueryResults$Citation[i] <- paste(Meta(gse_raw)$citation, collapse=" ")
  
  QueryResults$PMID[i] <- paste(Meta(gse_raw)$pubmed_id, collapse=" ")
  
  QueryResults$Contributor[i] <- paste(Meta(gse_raw)$contributor, collapse = " ")
  
  QueryResults$Date[i] <- paste(Meta(gse_raw)$submission_date, collapse= " ")
  
  QueryResults$Abstract[i] <- Meta(gse_raw)$summary
  
  rm(gse_raw)
}

#Getting an overview of our Query Result object with its new additions:

str(QueryResults)

#Adding empty columns to hold additional information that we will find while reviewing the dataset records:

QueryResults$Tissue<-character(length=nrow(QueryResults))
QueryResults$DevelopmentalStage<-character(length=nrow(QueryResults))
QueryResults$ManipulatedVariables<-character(length=nrow(QueryResults))
QueryResults$Notes<-character(length=nrow(QueryResults))

#Phrases that indicate the variable manipulated:
# Subjects were treated with ___ and rna-sequencing performed after testing behaviorally for...
# Subjects were divided into groups and one group experienced...
# Subjects received one of two interventions...

#Lets add some empty columns for taking inclusion/exclusion notes:

QueryResults$ManipulationUnrelatedToTopic<-character(length=nrow(QueryResults))
QueryResults$WrongTissue<-character(length=nrow(QueryResults))
QueryResults$NotBulkDissection_ParticularCellTypeOrSubRegion<-character(length=nrow(QueryResults))
QueryResults$IncorrectDevelopmentalStage<-character(length=nrow(QueryResults))
QueryResults$NotFullTranscriptome<-character(length=nrow(QueryResults))
QueryResults$MetadataIssues_MissingInfo_Retracted_Duplicated<-character(length=nrow(QueryResults))

QueryResults$Excluded<-character(length=nrow(QueryResults))
QueryResults$WhyExcluded<-character(length=nrow(QueryResults))

#Output the query results as a comma-separated variable file:

write.csv(QueryResults, "QueryResults.csv")

#You can see where you outputted the file by checking on the identity of your working directory.
getwd()
# [1] "/Users/hagenaue/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2026_TeamStress/R_Output_And_Results/Practice_Megan"


#########################
#########################

#How to extract the sample metadata in data frame format for a single GEO series:

#Use GEOQuery to pull down the full record (sample metadata and expression data)

library(GEOquery)

gse_raw <- getGEO("GSE237890", GSEMatrix=TRUE)


#You can see all of the goodies stashed in this object using the structure function:
str(gse_raw)

#Grab the "expression set" component of the object (item #1)
eset<-gse_raw[[1]]

#Grab the "phenoData" (sample metadata) for the expression set:
metadata_df <- pData(eset)

#You can see all of the components in the phenoData using str:
str(metadata_df)

#Or the column names for the phenoData:
colnames(metadata_df)

#Note: A lot of these column names have been auto named "characteristics_ch1..."

#To find out what those actually are, you can click on the data frame in your global environment (upper right)
#the actual variable name is stashed in the individual cells for the columns followed by a hyphen.

#Or you can view the first few rows of the phenoData using head:
head(metadata_df)

#Or you can write out the metadata as a .csv file and peruse it in a spreadsheet program:
getwd()

write.csv(metadata_df, "GSE237890_MetaData.csv")


#If we want to find out the distribution for our categorical variables (i.e., how many subjects are in each group), we can use tables and cross-tables:

#e.g.,
table(metadata_df$characteristics_ch1.3)

# social defeat stress: no social defeat stress (NIL) 
# 38 
# social defeat stress: social defeat stress (SD) 
# 42 

table(metadata_df$characteristics_ch1.4)

# adolescent environmental enrichment: social and environmental enrichment (EE) 
# 42 
# adolescent environmental enrichment: standard housing (NIL) 
# 38

#here's an example of a cross-table looking at how many subjects are in each subgroup defined by two variables (social defeat stress vs. no stress, environmental enrichment vs. standard housing):
table(metadata_df$characteristics_ch1.3, metadata_df$characteristics_ch1.4)

#here's an example of a cross-table looking at how many subjects are in each subgroup defined by three variables (social defeat stress vs. no stress, environmental enrichment vs. standard housing, hippocampus vs. nucleus accumbens):
table(metadata_df$characteristics_ch1.3, metadata_df$characteristics_ch1.4, metadata_df$source_name_ch1)



