#Exploring the distribution of gene expression values for the Gemma datasets:
#Megan Hagenauer
#7-10-2023

#This pulls up the log2 expression for the dataset (unfiltered):
GSE27532_Expression<-gemma.R::get_dataset_expression("GSE27532")
#Note that these values should probably range between -1 and 12 for RNA-Seq and 4-12 for microarray. 
#If you are only seeing values in a narrow range (e.g., between 2-3.58) that suggests there may be a problem, esp. if the data was external and no raw data was available
#I suspect this is happening because the external microarray data is often only available in log2 format and then accidentally log2-transformed again within Gemma's pipeline

str(GSE27532_Expression)

# Classes ‘data.table’ and 'data.frame':	25697 obs. of  20 variables:
# $ Probe                             : chr  "ILMN_1212607" "ILMN_1212612" "ILMN_1212619" "ILMN_1212628" ...
# $ GeneSymbol                        : chr  "Cradd" "Rcan2" "Mfap3l" "" ...
# $ GeneName                          : chr  "CASP2 and RIPK1 domain containing adaptor with death domain" "regulator of calcineurin 2" "microfibrillar-associated protein 3-like" "" ...
# $ NCBIid                            : chr  "12905" "53901" "71306" "" ...
# $ LA mouse saline hippocampus 4     : num  6.56 8.73 7.03 6.45 7.06 ...
# $ LA mouse saline hippocampus 3     : num  6.56 8.37 6.84 6.44 6.97 ...
# $ LA mouse saline hippocampus 2     : num  6.71 8.7 6.9 6.24 7.02 ...
# $ LA mouse saline hippocampus 1     : num  6.68 8.63 6.92 6.47 6.61 ...
# $ HA mouse saline hippocampus 4     : num  6.48 8.72 6.94 6.46 7.2 ...
# $ HA mouse saline hippocampus 3     : num  6.62 8.91 6.83 6.45 7.16 ...
# $ HA mouse saline hippocampus 2     : num  6.6 8.75 6.85 6.5 7.03 ...
# $ HA mouse saline hippocampus 1     : num  6.39 8.78 6.85 6.34 7.22 ...
# $ LA mouse desipramine hippocampus 4: num  6.62 8.49 6.8 6.45 7.06 ...
# $ LA mouse desipramine hippocampus 3: num  6.55 8.09 6.99 6.47 7.1 ...
# $ LA mouse desipramine hippocampus 2: num  6.61 8.33 6.83 6.5 7.01 ...
# $ LA mouse desipramine hippocampus 1: num  6.73 8.35 7.06 6.49 6.91 ...
# $ HA mouse desipramine hippocampus 4: num  6.51 8.58 6.86 6.43 7.05 ...
# $ HA mouse desipramine hippocampus 3: num  6.69 8.64 7.05 6.44 6.98 ...
# $ HA mouse desipramine hippocampus 2: num  6.51 8.21 6.94 6.47 7.23 ...
# $ HA mouse desipramine hippocampus 1: num  6.6 8.49 6.91 6.46 7.01 ...
# - attr(*, ".internal.selfref")=<externalptr> 
#   - attr(*, "call")= chr "https://gemma.msl.ubc.ca/rest/v2/datasets/GSE27532/data?filter=false"

#Pulling out the expression data and putting it in matrix form (which is a format that only includes 1 data type)
GSE27532_Expression_Matrix<-as.matrix(GSE27532_Expression[,-c(1:4)])

str(GSE27532_Expression_Matrix)
# num [1:25697, 1:16] 6.56 8.73 7.03 6.45 7.06 ...
# - attr(*, "dimnames")=List of 2
# ..$ : NULL
# ..$ : chr [1:16] "LA mouse saline hippocampus 4" "LA mouse saline hippocampus 3" "LA mouse saline hippocampus 2" "LA mouse saline hippocampus 1" ...

#Now I have a matrix of all numeric expression
#What are the range of values in this matrix?

summary(GSE27532_Expression_Matrix)

#This shows a histogram of all of the gene expression values in the dataset
#For Log2 expression values from microarray, these values will tend to range between 4-14
#In the original units (normalized fluorescent hybridization signal), that would be a minimum around 2^4:
2^4
#[1] 16
#... and a maximum around 2^14
2^14
#[1] 16384

#For Log2 expression values from RNA-Seq, these values will tend to range between 0-14.
#In non-log2 units, that would be a range of...
2^0
#[1] 1 
#Note - this is actually an artifact of analysis procedure - since 0 can't be log transformed, many protocols add 1 to all values before the transformation
2^14
#[1] 16384
#RNA-Seq units can be either counts per million (cpm), RPKM (Reads Per Kilobase Million), FPKM (Fragments Per Kilobase Million) or TPM (Transcripts Per Kilobase Million)
#These units control for different technical factors (total reads/sequencing depth, gene/transcript length)
https://www.reneshbedre.com/blog/expression_units.html
#CPM only corrects for sequencing depth/library size/total reads
#RPKM/FPKM corrects for transcript length and sequencing depth, but is primarily recommended for comparing the expression of different genes within a sample
#RPKM=single read RNA-Seq, FPKM=paired read RNA-Seq
#TPM also corrects for transcript length and sequencing depth
#Gemma uses... (I just e-mailed to ask - my guess is cpm)
  
#Here is a way to view a histogram of the log2 expression values for the entire sample: 
hist(GSE27532_Expression_Matrix)

#This shows an overview of the distribution of log2 expression values for each sample (column):
summary(GSE27532_Expression_Matrix)
# LA mouse saline hippocampus 4 LA mouse saline hippocampus 3 LA mouse saline hippocampus 2 LA mouse saline hippocampus 1
# Min.   : 6.192                Min.   : 6.192                Min.   : 6.192                Min.   : 6.192               
# 1st Qu.: 6.505                1st Qu.: 6.505                1st Qu.: 6.505                1st Qu.: 6.505               
# Median : 6.693                Median : 6.693                Median : 6.693                Median : 6.693               
# Mean   : 7.299                Mean   : 7.299                Mean   : 7.299                Mean   : 7.299               
# 3rd Qu.: 7.647                3rd Qu.: 7.647                3rd Qu.: 7.647                3rd Qu.: 7.647               
# Max.   :14.902                Max.   :14.902                Max.   :14.902                Max.   :14.902               
# HA mouse saline hippocampus 4 HA mouse saline hippocampus 3 HA mouse saline hippocampus 2 HA mouse saline hippocampus 1
# Min.   : 6.192                Min.   : 6.192                Min.   : 6.192                Min.   : 6.192               
# 1st Qu.: 6.505                1st Qu.: 6.505                1st Qu.: 6.505                1st Qu.: 6.505               
# Median : 6.693                Median : 6.693                Median : 6.693                Median : 6.693               
# Mean   : 7.299                Mean   : 7.299                Mean   : 7.299                Mean   : 7.299               
# 3rd Qu.: 7.647                3rd Qu.: 7.647                3rd Qu.: 7.647                3rd Qu.: 7.647               
# Max.   :14.902                Max.   :14.902                Max.   :14.902                Max.   :14.902               
# LA mouse desipramine hippocampus 4 LA mouse desipramine hippocampus 3 LA mouse desipramine hippocampus 2
# Min.   : 6.192                     Min.   : 6.192                     Min.   : 6.192                    
# 1st Qu.: 6.505                     1st Qu.: 6.505                     1st Qu.: 6.505                    
# Median : 6.693                     Median : 6.693                     Median : 6.693                    
# Mean   : 7.299                     Mean   : 7.299                     Mean   : 7.299                    
# 3rd Qu.: 7.647                     3rd Qu.: 7.647                     3rd Qu.: 7.647                    
# Max.   :14.902                     Max.   :14.902                     Max.   :14.902                    
# LA mouse desipramine hippocampus 1 HA mouse desipramine hippocampus 4 HA mouse desipramine hippocampus 3
# Min.   : 6.192                     Min.   : 6.192                     Min.   : 6.192                    
# 1st Qu.: 6.505                     1st Qu.: 6.505                     1st Qu.: 6.505                    
# Median : 6.693                     Median : 6.693                     Median : 6.693                    
# Mean   : 7.299                     Mean   : 7.299                     Mean   : 7.299                    
# 3rd Qu.: 7.647                     3rd Qu.: 7.647                     3rd Qu.: 7.647                    
# Max.   :14.902                     Max.   :14.902                     Max.   :14.902                    
# HA mouse desipramine hippocampus 2 HA mouse desipramine hippocampus 1
# Min.   : 6.192                     Min.   : 6.192                    
# 1st Qu.: 6.505                     1st Qu.: 6.505                    
# Median : 6.693                     Median : 6.693                    
# Mean   : 7.299                     Mean   : 7.299                    
# 3rd Qu.: 7.647                     3rd Qu.: 7.647                    
# Max.   :14.902                     Max.   :14.902 

#You'll note that they are all the same for this dataset. 
#That is because it is standard practice to quantile normalize microarray data:
https://en.wikipedia.org/wiki/Quantile_normalization
#Quantile normalization makes it so that all samples have the same overall distribution of values, but the original ranks of the datapoints are preserved for each sample
#Using quantile normalization reduces the influence of large-scale technical factors (e.g., total cDNA pippetted for each sample)
#Microarray data typically also is corrected for background fluorescence
#Together, the procedure of summarizing the data for the probes representing genes in a dataset, background correction, and quantile normalization is all part of a general pipeline
#This general pipeline is called RMA - robust multiarray analysis

GSE27532_Expression_SDperGene<-apply(GSE27532_Expression_Matrix, 1, sd)
hist(GSE27532_Expression_SDperGene)
#Mostly between 0-0.5, some values between 0.5-1.5
#Remember that a unit of 1 is a doubling.


#Let's try another one - this one is an external Illumina dataset: 
GSE73798_Expression<-gemma.R::get_dataset_expression("GSE73798")
GSE73798_Expression_Matrix<-as.matrix(GSE73798_Expression[,-c(1:4)])

hist(GSE73798_Expression_Matrix)
#Range 7-16
#Looks good - missing quite a few low values though. 
#I'm going to guess that low values were already filtered out in the external source

GSE73798_Expression_SDperGene<-apply(GSE73798_Expression_Matrix, 1, sd)
hist(GSE73798_Expression_SDperGene)
#Mostly between 0-0.5, some values as high as >1.5

#This one is RNA-Seq:
GSE81672_Expression<-gemma.R::get_dataset_expression("GSE81672")
GSE81672_Expression_Matrix<-as.matrix(GSE81672_Expression[,-c(1:4)])

hist(GSE81672_Expression_Matrix)
#Range -6-13
#So low values have definitely not been filtered out.
#Also means that the normalization method hasn't applied a simple +1 before log transformation
2^-6
#[1] 0.015625
#Since the values don't go lower than this, this must have been the value artificially added prior to log transformation.

#For any particular gene, here is the typical standard deviation:

GSE81672_Expression_SDperGene<-apply(GSE81672_Expression_Matrix, 1, sd)
hist(GSE81672_Expression_SDperGene)
#Mostly between 0-2, some values between 2-4
#Remember that a unit of 1 is a doubling.

#Here is an Agilent microarray dataset. 
#To analyze Agilent datasets it is necessary to use proprietary software
#So Gemma is dependent on the preprocesssed gene expression values provided "externally" by the original experimenters 
#We discovered last year that these datasets often look like they have been log2 transformed twice

GSE84183_Expression<-gemma.R::get_dataset_expression("GSE84183")
GSE84183_Expression_Matrix<-as.matrix(GSE84183_Expression[,-c(1:4)])

hist(GSE84183_Expression_Matrix)
#Values range from ~2.25-4.5
#In the "original units":
2^2.25
#[1] 4.756828
2^4.25
#[1] 19.02731
#... which looks a lot like a typical microarray range

#Looking at the SD:
GSE84183_Expression_SDperGene<-apply(GSE84183_Expression_Matrix, 1, sd)
hist(GSE84183_Expression_SDperGene)
#Mostly between 0-0.2, some values as high as 0.6


#That doesn't seem like quite enough to run on, so I went back to the original GEO record
#They call their gene expression data "normalized signal intensity"
#That is extremely ambiguous. It looks log transformed
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2228639

#So I looked at other records of data using the same platform. 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13912
#And it is totally the wild west.
#Folks are providing all sorts of crazy values as the normalized output to go with their experiments
#... and only some of them are labeled in detail regarding what they represent.

#But I did find two examples that were labeled Log2 signal output from two other experiments from the same platform 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM889978
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1121373
#As expected it looks like the values range between 1-16

#So I think the GSE84183 dataset was log2 transformed twice in the Gemma database

#If we want to reverse the log2 transformation:

GSE84183_Expression_Matrix_Corrected<-2^GSE84183_Expression_Matrix
hist(GSE84183_Expression_Matrix_Corrected)
#That looks more like typically microarray values
#Huge floor effect right below 5 - that must be unmeasurable with this experiment?

GSE84183_Expression_Corrected_SDperGene<-apply(GSE84183_Expression_Matrix_Corrected, 1, sd)
hist(GSE84183_Expression_Corrected_SDperGene)
#Now back to mostly ranging between 0-1, some values as high as 2.5

#Here's another Agilent dataset:
#On GEO it is already labeled Log2
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1501787

GSE61301_Expression<-gemma.R::get_dataset_expression("GSE61301")
GSE61301_Expression_Matrix<-as.matrix(GSE61301_Expression[,-c(1:4)])

hist(GSE61301_Expression_Matrix)
#These values range from -9 to 8.
#Interesting. There are many low values, but the values do reach almost more typical microarray levels.
#The distribution mostly looks shifted. Maybe there was some sort of background subtraction?

#Looking at the SD:
GSE61301_Expression_SDperGene<-apply(GSE61301_Expression_Matrix, 1, sd)
hist(GSE61301_Expression_SDperGene)
#Ranges mostly from 0-1, some values up to 4
#This dataset is probably actually o.k.
