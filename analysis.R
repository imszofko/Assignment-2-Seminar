##Install CRAN packages
#install.packages("openxlsx")

##Install Bioconductor packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install(c("edgeR", "clusterProfiler", "enrichplot", "org.Hs.eg.db"))

#Load the required libraries
library(edgeR)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(openxlsx)

#Import count data
counts <- read.table("E-MTAB-2523.counts.txt", header = TRUE, row.names = 1, sep = "\t")
print("Count data imported successfully.")
print(head(counts))

#Import sample data
sample_table <- read.table("E-MTAB-2523_sample table.txt", header = TRUE,  sep = "\t")
print("Sample data imported successfully.")
print(head(sample_table))

#Create a DGEList object
dge <- DGEList(counts = counts)
print("DGEList object created.")

#Calculate counts per million (CPM)
cpm <- cpm(dge)
print("CPM calculated.")

#Define a threshold for filtering
threshold <- 1
print(paste("Threshold for filtering set to", threshold))

#Keep genes with CPM > threshold in at least 2 samples
keep <- rowSums(cpm > threshold) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
print("Genes filtered.")

#Normalize the data
dge <- calcNormFactors(dge)
print("Data normalized.")

#Save filtered count data
write.table(dge$counts, file = "filtered_counts.txt", sep = "\t", quote = FALSE, col.names = NA)
print("Filtered count data saved.")