


pca <- read.table(paste0("covFiles/vsmcEuro_genoPCs.txt"), check.names=FALSE)

sex <- read.delim("genoFiles/vsmcEuro.fam", header = F)
sex$V3 <- paste0(sex$V1, "_", sex$V2)
sex <- data.frame(t(sex[, c(3,5)]))
colnames(sex) <- sex[1, ]
sex <- sex[-1, ]
sex <- sex[, which(names(sex) %in% names(pca))]
sex <- cbind(SampleID = "Sex", sex)

cov <- rbind(sex, pca)

write.table(cov, paste0("covFiles/vsmcEuro_Sex_10pc.txt"), sep = "\t", row.names = F, quote = F)

cat("Covariates file Ready  \n")
