#Example Generic Code for Making Simple Bivariate Plots & Quickly Getting Stats

#Replace this code with whatever variables you are actually interested in.
#These variables should either from the same object (e.g., the same SummarizedExperiment Object)
#Or they could be from different objects, but have the exact same samples in the same order
CategoricalVariable1<-colData(SummarizedExperiment_Filtered[[1]])$'organism part'
CategoricalVariable2<-colData(SummarizedExperiment_Filtered[[1]])$phenotype
NumericVariable1<-colData(SummarizedExperiment_Filtered[[1]])$LibrarySize
NumericVariable2<-TempExpressionData[1,]

#Examining the relationship between a numeric and categorical variable:
boxplot(NumericVariable1~CategoricalVariable1)
summary.lm(lm(NumericVariable1~CategoricalVariable1))

#Examining the relationship between a numeric and numeric variable:
plot(NumericVariable1~NumericVariable2)
summary.lm(lm(NumericVariable1~NumericVariable2))

#Examining the relationship between a categorical and categorical variable:
fisher.test(table(CategoricalVariable1, CategoricalVariable2))

