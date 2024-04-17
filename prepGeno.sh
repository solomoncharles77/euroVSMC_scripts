

module load plink2
module load bcftools
module load samtools
module load tabix



# Extract genotype data
plink2 --bfile ../genotypeData/rawGeno/allRawGeno/All_Imputed_Samples_All_Imputed_SNPs --keep ../genotypeData/rawGeno/eurGeno/Imputed_Data_TOPMed_hg38_European-only.fam --make-bed --out genoFiles/vsmcEuro --allow-extra-chr

# QC genotype data
plink2 --bfile genoFiles/vsmcEuro --geno 0.05 --hwe 1e-06 --maf 0.01 --recode vcf --out genoFiles/vsmcEuro --allow-extra-chr

bcftools view genoFiles/vsmcEuro.vcf --threads 6 -Oz -o genoFiles/vsmcEuro.vcf.gz
bcftools annotate -x INFO genoFiles/vsmcEuro.vcf.gz | bcftools +fill-tags  > genoFiles/vsmcEuro_QCed_tem.vcf

bcftools query -l genoFiles/vsmcEuro.vcf.gz > genoFiles/vsmcEuro_samples.txt

bcftools view genoFiles/vsmcEuro_QCed_tem.vcf -S genoFiles/vsmcEuro_samples.txt | sed 's/chr//g' | bgzip > genoFiles/vsmcEuro_QCed.vcf.gz
bcftools index -t genoFiles/vsmcEuro_QCed.vcf.gz
rm genoFiles/vsmcEuro_QCed_tem.vcf

# Prep gene expression bed file ---------------------------
Rscript euroVSMC_scripts/prepPheno.R

bgzip phenoFiles/geneExpr_qtlTools_Ready.bed
tabix -f phenoFiles/geneExpr_qtlTools_Ready.bed.gz

# Prep coveriates
plink2 --vcf genoFiles/vsmcEuro.vcf.gz --pca --out covFiles/vsmcEuro --allow-extra-chr 
Rscript euroVSMC_scripts/autoPrepGenoPCA.R
Rscript euroVSMC_scripts/prepCov.R


for chr in {1..22}
do
    echo "Processing chromosome ${chr}"
    # vcf files
    bcftools view genoFiles/vsmcEuro_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz
    
    # exp files
    zgrep -w "^#chr" phenoFiles/geneExpr_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/geneExpr_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed
    
    QTLtools cis --vcf genoFiles/chr${chr}.vcf.gz --bed phenoFiles/chr${chr}.bed.gz --cov covFiles/vsmcEuro_Sex_10pc.txt --nominal 1 --normal --std-err --out euroVSMC_cisEQTL/chr${chr}_euroVSMC_cisEQTL.txt.gz

    
done
