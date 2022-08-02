# Proteomics-and-Metabolomics-Integration
This is to test multi-omics integration tools with a proteomics and metabolomics example datasets.

The omics datasets used in this testing 'project' are also used in the article: Harutyunyan et al. 2022 Int. J. Mol. Sci. (https://doi.org/10.3390/ijms23116063).

Proteomics data is available via ProteomeXchange with dataset identifier PXD033987 

Metabolomics data is openly available at https://store.erc.monash.edu/experiment/view/15139/

Code used for testing will be added as the tools are tested out.

So far, MOFA appears to be the best tool to use for biomarker identification. It allows users to plot heatmaps and visualise the weight of each feature for individual samples. We can also extract data for downstream analysis.


It may also be useful to find methods that will allow us to identify the relationship between each feature per group. The more exploratory approaches such as MCIA and COMBI provide this option. However, it is unclear whether or not we can extract the values of these relationships.

To do:
1. find other approaches (maybe don't depend on a tool someone else wrote)
2. deep fusion/deep learning approaches could also be tried
