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
library(DESeq2)
#Import count data
countTable <- read.table("E-MTAB-2523.counts.txt", 
                     header = TRUE, 
                     row.names = 1, 
                     sep = "\t")
print("Count data imported successfully.")
print(head(countTable))

#Import sample data
sampleTable <- read.table("E-MTAB-2523_sample table.txt", 
                           header = TRUE,  
                           sep = "\t",
                           row.names = 1)
print("Sample data imported successfully.")
print(head(sampleTable))

#Create a factor for the carcinome and healthy patients

disease <- factor(sampleTable$disease)
sex <- factor(sampleTable$sex)
individual <- factor(sampleTable$individual)
'
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
'

'
I would approach the filtering differently, i would do a rowMeans(log2(cpm(countTable) + 1)) for filtering and then save the 
filtered genes that are greater than one counts = counts[meanlog2CPM > 1, ]
This is how I remember it and usually do it and it is typically no issue

As for the normalization, I do not think we need to normalize it before creating the dds thing using the Deseqmatrix fucntion. 
it can be normalized with the deseq function after
'

#filtering out the low count genes in the dataset
meanlog2CPM <- rowMeans(log2(cpm(countTable + 1)))

#removing values that are lower than 1 by keeping row values that are higher than 1
countTable <- countTable[meanlog2CPM > 1, ]

##Perform statistical analysis to identify the DEGs between the carcinoma and the normal samples
##Genes should be filtered on FDR and log2FC

dds <- DESeqDataSetFromMatrix(countData = countTable,
                             colData = sampleTable,
                             design = ~disease)
#Compute the statistical analysis by the function DESeq
dds2 <- DESeq(dds)

#Results for dds2
results <- results(dds2)
summary(results) #Shows DEGs, low and high regulated, and how many outliers there are

#Statistical analysis part 
'Looking back at previous things ive worked on, I think that this part fits best here'

#Create a DGEList object
dge <- DGEList(counts = countTable)
print("DGEList object created.")

#Normalize the data
dge <- calcNormFactors(dge)
print("Data normalized.")

studyDesign = data.frame(Individual = individual, Disease = disease)
desMatrix = model.matrix(~ + Disease, data = studyDesign)
colnames(desMatrix)[1:2] = c('Normal','Carcinoma')

#Fit linear model, using the normalized data and designmatrix
dge <- estimateDisp(dge, desMatrix, robust = T)
fit <- glmQLFit(dge, desMatrix, robust = T)

contrastMatrix = makeContrasts(CarcinomaVsNormal = Carcinoma - Normal,
                               levels = desMatrix)
                               
#Stastical testing
qlfRes <- glmQLFTest(fit, contrast = contrastMatrix)
topRes <- topTags(qlfRes, n = nrow(fit$counts))
topRes <- subset(topRes$table, abs(logFC) > 2 & FDR < 0.05)

#Just to make sure they are un order
sorted = topRes[order(topRes$logFC, decreasing = T),]
dim(topRes) #12979 genes
print(topRes[1:10,])
topRes

#Save filtered count data
write.table(dge$counts, file = "filtered_counts.txt", sep = "\t", quote = FALSE, col.names = NA)
print("Filtered count data saved.")

#Save genes filtered in FDR and Log2FC
write.table(topRes, file = "DEGs.txt", sep = '\t',
            quote = F, col.names = NA)

##Perform ORA-analysis with the functions enrichGO and enrichKEGG in the clutserProfile package
##



