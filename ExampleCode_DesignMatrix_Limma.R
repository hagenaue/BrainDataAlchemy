#Making a design matrix:

table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
# imipramine                        ketamine reference substance role,saline 
# 7                               6                               3
table(SummarizedExperiment_Subset_noBad_Filtered$phenotype)
#Non-responder,susceptible toward               susceptible toward     susceptible toward,Responder 
#7                                3                                6 

str(SummarizedExperiment_Subset_noBad_Filtered$treatment)
#chr [1:16]
#currently a character vector

SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-as.factor(SummarizedExperiment_Subset_noBad_Filtered$treatment)
levels(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
#[1] "imipramine"                      "ketamine"                        "reference substance role,saline"

SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-relevel(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor, ref="reference substance role,saline")
levels(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
#[1] "reference substance role,saline" "imipramine"                      "ketamine" 

str(SummarizedExperiment_Subset_noBad_Filtered$phenotype)
#chr [1:16]
#currently a character vector

library(limma)

colData(SummarizedExperiment_Subset_noBad_Filtered)$treatment_factor

design <- model.matrix(~treatment_factor+LibrarySize, data=colData(SummarizedExperiment_Subset_noBad_Filtered))
design

#testing it out with the data for 1 gene 
#we have to tell the lm function to not add an intercept
summary.lm(lm(ExpressionData_Subset_noBad_Filtered[1,]~0+design))


#Applying the model to all genes using limma:

fit <- lmFit(ExpressionData_Subset_noBad_Filtered, design)

#Adding an eBayes correction to help reduce the influence of outliers/small sample size on estimates
efit <- eBayes(fit, trend=TRUE)

dt<-decideTests(efit)
summary(decideTests(efit))
# (Intercept) treatment_factorimipramine treatment_factorketamine LibrarySize
# Down          4320                          0                       59           0
# NotSig        5135                      21686                    21512       21686
# Up           12231                          0                      115           0

str(efit)

#Writing out the results into your working directory:
write.fit(efit, adjust="BH", file="Limma_results.txt")
write.csv(rowData(SummarizedExperiment_Subset_noBad_Filtered), "Annotation_LimmaResults.csv")
