#Importing the count data
countTable <- read.table(file = "data/counts.tsv",
                         header = TRUE,
                         as.is = TRUE,
                         sep = "\t",
                         row.names = 1)
# Check the dimensions (number of rows (genes) and columns (samples)) of countTable using dim()
dim(countTable)


# Define experimental groups (Dose levels) and associate them with sample names
Dose <- factor(c(rep("Low", 5), rep("High", 10)), levels = c("Low", "High"))
sample_names <- colnames(countTable)
stopifnot(identical(sample_names, c(paste0("Low_", 1:5), paste0("High_", 1:10))))
#Filter low counts
library(edgeR)
library(limma)
meanLog2CPM <- rowMeans(log2(cpm(countTable) + 1))
# Visualize distribution of meanLog2CPM by plotting a histogram
hist(meanLog2CPM)
# Extract genes with meanLog2CPM > 1
countTable <- countTable[meanLog2CPM > 1, ]
dim(countTable)


#Define DGEList
dge <- DGEList(countTable)
dge$samples

dge <- calcNormFactors(dge)
dge$samples


#Define design matrix
designDF <- data.frame(Dose = Dose)
design <- model.matrix(~ 0 + Dose , data = designDF)


#Estimate dispersion parameters
dge <- estimateDisp(dge, design, robust = TRUE)


#Fit model
fit <- glmQLFit(dge, design, robust = TRUE)


# Inspect the coefficients for the first six genes using head() function
head(fit$coefficients)

#Define contrast matrix
contrastMat <- makeContrasts(HighVsLow = DoseHigh - DoseLow,levels = design)
contrastMat

#Perform hypothesis testing
qlfRes <- glmQLFTest(fit, contrast = contrastMat)
#Statistical testing
topRes <- topTags(qlfRes, n = nrow(fit$counts))


#Identifing differentially expressed genes
topRes <- subset(topRes$table, abs(logFC) > 1 & FDR < 0.05)
nrow(topRes)

# Order by FDR (ascending) and take the top 10 genes
top10 <- head(topRes, 10)

# Print the top 10 genes in a nice table
print(top10)

write.csv(topRes, "results/differential_expression.csv")
write.csv(top10, "results/top_genes.csv")
