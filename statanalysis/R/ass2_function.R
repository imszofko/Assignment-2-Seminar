#' Assignment 2 Statistical Analysis Function
#' 
#' The following package will input an RNA-seq count data and sample data into R.
#'   
#' @param countTable,sampleData The name/path of the file in quotes marks
#' @return Count data table
#' @examples
#' count_data <- countData("E-MTAB-2523.counts.txt");
#' 
#' @export
countData <- function(countTable){
  # Load necessary libraries
  library(edgeR)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(openxlsx)
  library(DESeq2)
  library(AnnotationDbi)
  library(ggplot2)
  
  # Import count data
  countTable <- read.table("E-MTAB-2523.counts.txt", header = TRUE, row.names = 1, sep = "\t")
  print("Count data imported successfully.")
  print(head(countTable))
}

#' Sample data reading
#' 
#' Reads in the data of the sample table and creates factors for the information
#' @param sampleData The name/path to the sample data file
#' @return The sample data and the factors
#' @examples
#' sample_data <- sampleData("E-MTAB-2523_sample table.txt");
#' @export
sampleData <- function(sampleData){
  # Import sample data
  sampleTable <- read.table("E-MTAB-2523_sample table.txt", header = TRUE, sep = "\t", row.names = 1)
  print("Sample data imported successfully.")
  print(head(sampleTable))
  
  # Create factors for the carcinoma and healthy patients
  disease <- factor(sampleTable$disease)
  sex <- factor(sampleTable$sex)
  individual <- factor(sampleTable$individual)
  return(sampleTable, disease, sex, individual);
}

#' Performing filtering
#' 
#' Filters by rowmeans(log2(cpm())+1)
#' @param countTable Filter count data
#' @return Count table with low filtered values removed
#' @examples
#' filterdata <- filtering(countTable);
#' @export
filtering <- function(countTable){
  # Filter out low count genes
  meanlog2CPM <- rowMeans(log2(cpm(countTable + 1)))
  countTable <- countTable[meanlog2CPM > 1, ]
  return(countTable);
}

#' Statistical Analysis
#'
#' Runs Statistical Analysis on the count data
#' @param countTable,sampleTable Variable names for the count data and sample table
#' @return dge and desMatrix
#' @examples
#' stats <- statanalysis(countTable, sample_data);
#' @export 
statanalysis <- function(countTable, sampleTable){
  # Perform statistical analysis to identify DEGs
  dds <- DESeqDataSetFromMatrix(countData = countTable, colData = sampleTable, design = ~disease)
  dds2 <- DESeq(dds)
  results <- results(dds2)
  summary(results)
}
#' DGE function
#' 
#' Create DGE table
#' @param countTable countTable after filtering
#' @return dge table
#' @examples
#' dge <- dge(countTable);
#' @export
dgeObj <- function(countTable){
  
  # Create a DGEList object and normalize the data
  dge <- DGEList(counts = countTable)
  dge <- calcNormFactors(dge)
  print("Data normalized.")
}

#' Linear Models
#' 
#' Makes design matrix and Fits linear models and performs stat testing
#' @param dge dge variable and variable for the design matrix 
#' @return Data files and png files for enrichment analysis visual representation
#' @examples
#' linearmodel <- linmodel(dge);
#' @export
linmodel <- function(dge, desMatrix){
  
  studyDesign <- data.frame(Individual = individual, Disease = disease)
  desMatrix <- model.matrix(~ Disease, data = studyDesign)
  colnames(desMatrix)[1:2] <- c('Normal', 'Carcinoma')
  
  # Fit linear model and perform statistical testing
  dge <- estimateDisp(dge, desMatrix, robust = TRUE)
  fit <- glmQLFit(dge, desMatrix, robust = TRUE)
  contrastMatrix <- makeContrasts(CarcinomaVsNormal = Carcinoma - Normal, levels = desMatrix)
  qlfRes <- glmQLFTest(fit, contrast = contrastMatrix)
  topRes <- topTags(qlfRes, n = nrow(fit$counts))
  topRes <- subset(topRes$table, abs(logFC) > 2 & FDR < 0.05)
  sorted <- topRes[order(topRes$logFC, decreasing = TRUE),]
  print(topRes[1:10,])
  
  # Save filtered count data and DEGs
  write.table(dge$counts, file = "filtered_counts.txt", sep = "\t", quote = FALSE, col.names = NA)
  print("Filtered count data saved.")
  write.table(topRes, file = "DEGs.txt", sep = '\t', quote = FALSE, col.names = NA)
  print("DEGs saved.")
  
  # Perform ORA-analysis with enrichGO and enrichKEGG
  changeNames <- mapIds(org.Hs.eg.db, keys = rownames(topRes), keytype = "SYMBOL", column = "ENTREZID", multiVals = 'first')
  topRes$EntrezIds <- changeNames
  duplicated <- topRes$EntrezIds[duplicated(topRes$EntrezIds)]
  paste0('Number of NA and duplicates: ', length(duplicated))
  duplicatedGeneEntrez <- duplicated[!is.na(duplicated)]
  paste0('Number of Gene Entrez Ids: ', length(duplicatedGeneEntrez))
  NAs <- duplicated[is.na(duplicated)]
  paste0('Number of NAs: ', length(NAs))
  geneEntrez <- topRes$EntrezIds[!is.na(topRes$EntrezIds)]
  
  # Perform ORA-analysis with enrichGO
  print("Performing ORA-analysis with enrichGO...")
  ego <- enrichGO(gene = topRes$EntrezIds, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
  print("enrichGO analysis complete.")
  print(head(ego))
  
  # Perform ORA-analysis with enrichKEGG
  print("Performing ORA-analysis with enrichKEGG...")
  ekegg <- enrichKEGG(gene = topRes$EntrezIds, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)
  print("enrichKEGG analysis complete.")
  print(head(ekegg))
  
  # Export DEGs to Excel
  print("Exporting DEGs to Excel...")
  write.xlsx(topRes, file = "DEGs.xlsx")
  print("DEGs exported to Excel.")
  
  # Create plots to visualize the results of the enrichment analysis
  print("Creating cnetplot for enrichGO results...")
  cnetplot(ego, showCategory = 10)
  ggsave("cnetplot_enrichGO.png")
  
  print("Creating dotplot for enrichGO results...")
  dotplot(ego, showCategory = 10)
  ggsave("dotplot_enrichGO.png")
  
  print("Creating treeplot for enrichGO results...")
  ego_sim <- pairwise_termsim(ego)
  treeplot(ego_sim, showCategory = 10)
  ggsave("treeplot_enrichGO.png")
  
  print("Creating cnetplot for enrichKEGG results...")
  cnetplot(ekegg, showCategory = 10)
  ggsave("cnetplot_enrichKEGG.png")
  
  print("Creating dotplot for enrichKEGG results...")
  dotplot(ekegg, showCategory = 10)
  ggsave("dotplot_enrichKEGG.png")
  
  print("Creating treeplot for enrichKEGG results...")
  ekegg_sim <- pairwise_termsim(ekegg)
  treeplot(ekegg_sim, showCategory = 10)
  ggsave("treeplot_enrichKEGG.png")
}