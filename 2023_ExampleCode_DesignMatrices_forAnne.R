#Code for Anne for making design matrices -
#Two examples - Will need to be adapted to your subsetted filtered objects

################

#GSE172133

table(SummarizedExperiment[[1]]$treatment)
table(SummarizedExperiment[[1]]$timepoint)

SummarizedExperiment[[1]]$treatment_factor<-as.factor(SummarizedExperiment[[1]]$treatment)
levels(SummarizedExperiment[[1]]$treatment_factor)
#[1] "chronic constriction injury (CCI),nervous system injury" "reference subject role"  
SummarizedExperiment[[1]]$treatment_factor<-relevel(SummarizedExperiment[[1]]$treatment_factor, ref="reference subject role")
levels(SummarizedExperiment[[1]]$treatment_factor)
#[1] "reference subject role"                                  "chronic constriction injury (CCI),nervous system injury"

SummarizedExperiment[[1]]$timepoint_factor<-as.factor(SummarizedExperiment[[1]]$timepoint)
levels(SummarizedExperiment[[1]]$timepoint_factor)
#[1] "1 d"            "14 d"           "3 d"            "7 d"            "Not Applicable"
#This is going to be a co-variate in the model, with no interaction term, so the reference doesn't really matter

library(limma)

design <- model.matrix(~SummarizedExperiment[[1]]$treatment_factor+SummarizedExperiment[[1]]$timepoint_factor)

design
#dummy variable
#reference subject role=0
#injury=1
#Slope in the regression model is actually the group difference.

#timepoint
#1d= 0
#3d= 1
#can't set the other levels as 2, 3, 4, 5... because R will think it is a continuous variable (ranked)
#So we break it into multiple variables

#variable 1:
#1d=0
#3d=1
#everything else... also 0

#variable 2:
#1d=0
#7d=1
#everything else... also 0

#variable 3:
#1d=0
#14d=1
#everything else... also 0

#(Type III) Regression gives us the results for 1 variable while controlling for all others
#So we end up with each of these variables giving us just the level set at 1 vs. the reference in the output

#the design matrix:
design

#The Intercept is defined by the contrast structure
#Contrast treatment means that the intercept is at the reference level for all variables
#This is different from the default for ANOVA post-hoc contrast output

#This may crash due to timepoint not applicable being the same as no injury

colnames(design)
cbind(design[,2], design[,6])
table(design[,2], design[,6])
#    0  1
# 0  0  3
# 1 15  0

#demonstration with simple linear regression:
ExpressionData_Temp<-assay(SummarizedExperiment[[1]])
head(ExpressionData_Temp)

summary.lm(lm(ExpressionData_Temp[1,]~design))
#Seems unhappy that time point not applicable is same as injury

summary.lm(lm(ExpressionData_Temp[1,]~0+design[,-6]))

#get rid of that not applicable column that is causing issues
design<-design[,-6]

################

#GSE92718

rm(SummarizedExperiment)
SummarizedExperiment<-gemma.R::get_dataset_object("GSE92718", type = 'se', filter=TRUE, consolidate="average")
SummarizedExperiment
#dim: 16863 24 

table(SummarizedExperiment[[1]]$treatment, SummarizedExperiment[[1]]$timepoint) 
#                           2 w 8 w
# cuff operation             6   6
# reference substance role   6   6

SummarizedExperiment[[1]]$treatment_factor<-as.factor(SummarizedExperiment[[1]]$treatment)
levels(SummarizedExperiment[[1]]$treatment_factor)
#[1] "cuff operation"           "reference substance role" 
SummarizedExperiment[[1]]$treatment_factor<-relevel(SummarizedExperiment[[1]]$treatment_factor, ref="reference substance role")
levels(SummarizedExperiment[[1]]$treatment_factor)
#[1] "reference substance role" "cuff operation" 

SummarizedExperiment[[1]]$timepoint_factor<-as.factor(SummarizedExperiment[[1]]$timepoint)

design <- model.matrix(~SummarizedExperiment[[1]]$treatment_factor+SummarizedExperiment[[1]]$timepoint_factor)

design
