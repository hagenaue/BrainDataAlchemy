# BrainDataAlchemy
Code to mine for the summer Brain Data Alchemy meta-analysis projects.

## Pipeline from Summer 2022:

Now overviewed in detail in the published protocol:
Hagenauer M, Rhoads C, Xiong J, Nguyen DM, Hernandez E, Saffron A, Kondur A, Flandreau EI: Protocol: Brain Data Alchemy Project: Meta-Analysis of Re-Analyzed Public Transcriptional Profiling Data in the Gemma Database, 08/2024., dx.doi.org/10.17504/protocols.io.j8nlk84jxl5r/v1

The published protocol provides links to each of the relevant code documents in this folder.

## Pipeline from Summer 2023:

This code is all within the MetaAnalysis_GemmaDatasets folder:
https://github.com/hagenaue/BrainDataAlchemy/tree/main/MetaAnalysis_GemmaDatasets%20


### Gemma API: 

1) This is example code for using the Gemma API to conduct a systematic search for potential datasets to include in the meta-analysis

https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_Antidepressants_New_Code_forGemmaSearch.R


### Analyzing a single dataset:

2) This is example code for exploring and re-analyzing a single dataset:

https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_Example_ExploringADataSet_MoreDetail.R

3) For RNA-Seq datasets, we later realized that we would also need library size information - this code is a patch that adds that information right after the dataset is read into R:
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_Example_GettingLibrarySizeInfo.R

4) This is a simpler version of the example code for exploring and re-analyzing a single dataset (meant for microarray experiments):
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_Example_ExploringADataSet_SimplerForMicroarray.R

5) This is example code for navigating a situation where there isn’t a summarized experiment object available for a dataset:
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_DealingWithNoSummarizedExperiment_ForChristabel.R

6) Writing out Gemma summarized experiment objects in a format that is easily navigable by people unfamiliar with R:
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_WritingOutGemmaDatasets_forDocumentation.R


### Running a meta-analysis:

7) Reading in and formatting the limma results for the meta-analysis:

https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_Formatting_LimmaResults_forMetaAnalysis.R

8) Aligning results from different datasets and different species:
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_AligningLimmaResults_forMetaAnalysis.R

Note: That alignment code uses gene homolog/ortholog information that comes from these two files:
https://github.com/hagenaue/BrainDataAlchemy/blob/main/HOM_MouseVsRat.csv
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MouseVsRat_NCBI_Entrez_JacksonLab_20230727.csv

… and those files were downloaded from Jackson Labs (07/2023) and cleaned up using this code:
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_FormattingRatMouseOrthologDatabase.R

9) Running a simple (intercept only) version of the meta-analysis:
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_MetaAnalysisCode.R

10) Running a fancier version of the meta-analysis (intercept and modifying variable):
https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_ExampleCode_FancierMetaAnalysis_Christabel.R





