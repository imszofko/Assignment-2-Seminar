# Assignment-2-Seminar
Group Work Repo for the class Methods and Tools in Software Development for Bioinformatics

## First Objective 
Task in the mini-project is to create an R package that can perform enrichment analysis on an RNA-seq dataset. The package should be based on the R packages edgeR, clusterProfiler, enrichplot, org.Hs.eg.db and openxlsx.
One member in the group starts by creating a GitHub repo, committing a README file and adding the other group members as collaborators. In order to make use of GitHub’s pull request feature (one of the grading criteria), make sure to work in separate Git branches asnew features are added, and later merge these into the main branch. 
- The main goal for this first part is to understand how to create and set up a repository and afterwards to learn how to work together, on Github, with others when developing software.

## Second Objective
Objective is to create an R package with the group. The R package should be able to:
- Import RNA-seq count data into R and perform filtering to remove low-expressed genes. You will use the E-MTAB-2523 RNA-seq counts as a case study dataset when developing the package. This can be downloaded together with its sample table from the seminar page in Canvas. The dataset contains partially paired samples from colorectal carcinoma and health control tissue.


- Perform statistical analysis to identify differentially expressed genes between the carcinoma and normal samples. Genes should be filtered on FDR and log2 fold change.


- Perform over-representation analysis with the functions enrichGO and enrichKEGG in the clusterProfiler package. The key type for the genes in the dataset is SYMBOL. Note that enrichKEGG expects Entrez gene IDs (ENTREZID), so IDs need to be converted before analyzing KEGG pathways. The clusterProfiler package provides the convenient bitr function for this. As an alternative, you may perform gene set enrichment analysis with the gseGO and gseKEGG functions.


- Export the table with differentially expressed genes to Excel and create plots to visualize the results of the enrichment analysis. Examples of functions to use are cneplot, dotplot and treeplot. The last one depends on a term similarity matrix, which is calculated by calling pairwise_termsim on the enrichment analysis results.


### Contributors
Sophia Arias (email: b22sopar@student.his.se)

Natalia Medrano Do Nascimento (email: a22natme@student.his.se)
