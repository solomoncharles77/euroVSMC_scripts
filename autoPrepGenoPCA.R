

# read in the eigenvectors, produced in PLINK
eigenvec <- read.table("covFiles/vsmcEuro.eigenvec", header = FALSE, skip=0)
rownames(eigenvec) <- eigenvec[,1]
eigenvec <- eigenvec[,2:ncol(eigenvec)]
colnames(eigenvec) <- paste0("genoPC", c(1:ncol(eigenvec)))
eigenvec <- data.frame(t(eigenvec))
eigenvec <- data.frame(SampleID = rownames(eigenvec), eigenvec)
colnames(eigenvec) <- sub("X", "", colnames(eigenvec))
write.table(eigenvec, "covFiles/vsmcEuro_genoPCs.txt", sep = "\t", quote = F, col.names = colnames(eigenvec))

