
# # This script takes one argument -  the prefix name of the plink processed output
# aha <- commandArgs(trailingOnly = TRUE)
# 
# aha2 <- paste0("genoFiles/", aha, ".Geno.txt")
# aha3 <- paste0("phenoFiles/", aha, ".txt")

# load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))


# Import and process geno data --------------------------------------------
genoSamps <- read.table("genoFiles/vsmcEuro_samples.txt")
genoSamps$V2 <- sapply(genoSamps$V1, function(x) {unlist(strsplit(x, split = "_"))[2]})
head(genoSamps)

cat("\n")

cat("This is a view of the expression data \n")
cat("\n")
# Import and process expr data --------------------------------------------
expr <- data.frame(fread("/home/c/cs806/VSMC_Raw_gene_counts_matrix.txt.gz"))
colnames(expr) <- gsub("A", "S", colnames(expr))
# Harmonize names and samples in geno and expr ----------------------------


# Extract and normalize expr samples present in geno data
expr1 <- expr[, which(colnames(expr) %in% genoSamps$V2)]
rownames(expr1) <- expr$GeneID
expr1[1:10, 1:10]


cat("\n")
cat("Normalizing expression data \n")
cat("\n")

# normalize expr
sampleTable <- data.frame(condition = factor(rep(c("HUVEC"), ncol(expr1))))
rownames(sampleTable) <- colnames(expr1)

dds <- DESeqDataSetFromMatrix(round(expr1), sampleTable, ~1)
ddsF <- dds[ rowSums(counts(dds)) > ncol(dds), ]
vst=varianceStabilizingTransformation(ddsF)
normExpr <- as.data.frame(assay(vst))

# Rename to match name in vcf file, order samples and Re-add geneID columns
names(normExpr) <- genoSamps$V1[match(names(normExpr), genoSamps$V2)]
normExpr <- normExpr[, order(names(normExpr))]
normExpr <- cbind(geneID = rownames(normExpr), normExpr)

cat("This is a view of the normalized expression data \n")
cat("\n")
normExpr[1:10, 1:10]


# Prep bed for normalized genes -----------------------------------
# Add coordinate information 
exprCoord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords_withStrand.txt"))
colnames(exprCoord)[1] <- "geneID"
normExprCoord <- exprCoord[which(exprCoord$geneID %in% rownames(normExpr)), ]
normExprCoord <- merge(normExprCoord, normExpr, by = "geneID")
normExprCoord <- normExprCoord[, c(2:4,1,6,5, 7:ncol(normExprCoord))]
colnames(normExprCoord)[1:6] <- c("#chr", "start", "end", "gene", "length", "strand")
normExprCoord <- normExprCoord[order(normExprCoord[,1], normExprCoord[,2]), ]
normExprCoord[1:10, 1:10]
cat("\n")
cat("Exporting QTL ready files \n")
cat("\n")

# Export files for QTLtools ---------------------------------------------
exprFileName <- paste0("phenoFiles/geneExpr_qtlTools_Ready.bed")
fwrite(normExprCoord, file = exprFileName, sep = "\t")

# system(paste0("module load samtools"))
# system(paste0("module load tabix"))
# system(paste0("bgzip ", exprFileName))
# system(paste0("tabix ", exprFileName, ".gz"))
