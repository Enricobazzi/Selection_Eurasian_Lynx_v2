---
title: "RDA_pcs"
author: "Enrico"
date: "2/10/2022"
output: html_document
editor_options:
  chunk_output_type: console
---

### RDA-1 Preparing the data

In preparation for RDA, I generated a table of values for each sample of the selected principal components of environmental predictors (see 2-Preparing_Environmental_Data -> PCA of samples based on variables). I have generated two different tables, one with only the first two PCs and one with the first five. A decision between analyzing only 2 or 5 PCs will have to be made. For now I will analyze with a multivariate RDA both sets and with univariate BayPass each PC separately.

Samples from eroded populations (Balkans, Carpathians, Poland and Norway) as well as c_ll_vl_0137, c_ll_og_0184, c_ll_og_0187, c_ll_tu_0154 were excluded from selection scans.

To upload it to CESGA or genomics server:
```{bash}
# 2 pcs
scp ~/Documents/Selection_Eurasian_Lynx_v2/2-Prepare_Environmental_Data/twopcs_data_matrix.tsv \
csebdeba@ft2.cesga.es:/mnt/lustre/scratch/home/csic/ebd/eba/Selection_Eurasian_Lynx/RDA
# 5 pcs
scp ~/Documents/Selection_Eurasian_Lynx_v2/2-Prepare_Environmental_Data/fivepcs_data_matrix.tsv \
csebdeba@ft2.cesga.es:/mnt/lustre/scratch/home/csic/ebd/eba/Selection_Eurasian_Lynx/RDA

# 2 pcs
scp ~/Documents/Selection_Eurasian_Lynx_v2/2-Prepare_Environmental_Data/twopcs_data_matrix.tsv \
ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/RDA
# 5 pcs
scp ~/Documents/Selection_Eurasian_Lynx_v2/2-Prepare_Environmental_Data/fivepcs_data_matrix.tsv \
ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/RDA

```

I also need to convert my VCF file to RAW format, for the RDA function to interpret it. I will do this using PLINK v1.9, eliminating missing data (geno 0) and the unwanted samples
```{bash}
cd /mnt/lustre/scratch/home/csic/ebd/eba/Selection_Eurasian_Lynx/RDA

compute

# load plink module
module load cesga/2020
module load plink

# convert VCF of the Selection dataset (no bottleneck populations and maf>0.05) to RAW
plink --vcf ${STORE2}/lynx_genome/lynx_data/LyCaRef_vcfs/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--geno 0 \
--recode A --out ll_wholegenome_filtered_snps_finalset
```

### RDA-2 Running RDA

To run RDA I will open an R session in a node with high memory or use most free genomics server
```{bash}
# on cesga
cd /mnt/lustre/scratch/home/csic/ebd/eba/Selection_Eurasian_Lynx/RDA
compute --mem 40
module load cesga/2020
# on genomics
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/RDA

R
```

Once in R load necessary libraries (install missing ones first with install.packages("name_of_package")) and the function I will use to define outliers based on Standard Deviation
```{R}
# load libraries
library(adegenet)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(parallel)

# load function to define outliers based on SD
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]              # locus names in these tails
}
```
Read input files
```{R}
# load environmental data - 2 PCs
env.predictors <- read_tsv("twopcs_data_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")

# load environmental data - 5 PCs
env.predictors <- read_tsv("fivepcs_data_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")


# load SNPs data as genlight object
gt_data <- read.PLINK("ll_wholegenome_filtered_snps_finalset.raw")

# convert SNPs data to dataframe
gt_data_tsv <- data.frame(as.matrix(gt_data))
```

Run RDA
```{R}
rda <- rda(gt_data ~ ., data=env.predictors, scale=T)
```
To explore which axes of the RDA results were significant and will be further analyzed I first drew a screeplot to visualize inertia of each axis
```{R}
# create screeplot of RDA results - 2 PCs
pdf(file = paste0("screeplot_2pcs.pdf"), width = 4, height = 4)
screeplot(rda)
dev.off()

# create screeplot of RDA results - 5 PCs
pdf(file = paste0("screeplot_5pcs.pdf"), width = 4, height = 4)
screeplot(rda)
dev.off()
```
To check for significance I ran an anova of the canoncical correspondence analysis. For memory constraints only 100 permutations were investigated. This was run on the genomics EBD server, as RAM usage was too high for CESGA (300+ GB used)
```{R}
# calculate significance of RDA axes
options(mc.cores=1)
signif.axis <- anova.cca(rda, by="axis", parallel=getOption("mc.cores"), permutations = how(nperm=100), cutoff=0.05)
```
The following results were found:
```{R}
signif.axis
# Permutation test for rda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 100
# 
# Model: rda(formula = gt_data ~ PC1 + PC2 + PC3 + PC4 + PC5, data = env.predictors, scale = T)
#          Df Variance      F   Pr(>F)   
# RDA1      1   128181 5.3548 0.009901 **
# RDA2      1    85217 3.5600 0.009901 **
# RDA3      1    37721 1.5758 0.019802 * 
# RDA4      1    27614 1.1536 0.306931   
# RDA5      1    22114 0.9238            
# Residual 63  1508070                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
Significance value of third axis is low and very close to significance threshold, so we will consider it while looking for candidate SNPs for selection. MANUALLY INPUT ACCORDINGLY
```{R}
# number of significant axes - 2 PCs
nsig.axes <- 2
# number of significant axes - 5 PCs
nsig.axes <- 3
```
To investigate candidate SNPs I will first get the loadings of each SNP for each significant axis of the RDA and write a table of them for safekeeping :)
```{R}
# load all loadings of significant axes
load.rda <- scores(rda, choices=c(1:nsig.axes), display="species")

# write a table with all loadings of significant axes - 2 PCs
write.table(x = rownames_to_column(data.frame(load.rda), var="SNP"),
            file = "rda_2pcs_loadings_sigaxes.tsv",
            quote=FALSE,  col.names = T, row.names = F, sep= "\t")

# write a table with all loadings of significant axes - 5 PCs
write.table(x = rownames_to_column(data.frame(load.rda), var="SNP"),
            file = "rda_5pcs_loadings_sigaxes.tsv",
            quote=FALSE,  col.names = T, row.names = F, sep= "\t")
```
SNPs with loading values that exceeded 2.5 times the SD were chosen as candidates. This threshold is relatively low which means more true positives but also more false positives. As this list of candidates will be then intersected with other selection tests (BayPass), I felt more confident in being a bit more "forgiving" in this stage.
```{R}
# extract candidates of each axis
cand <- data.frame()
for(i in 1:nsig.axes){
  candN <- outliers(load.rda[,i],2.5)
  candN <- cbind.data.frame(rep(i,times=length(candN)), names(candN), unname(candN))
  colnames(candN) <- c("axis","snp","loading")
  cand <- rbind(cand, candN)
}
cand$snp <- as.character(cand$snp)

# remove duplicate entries (snps that are outliers for more than one axis)
cand <- cand[!duplicated(cand$snp),]
ncand <- NROW(cand)
```
In order to visually assign a predictor to each candidate SNP, we can calculate correlation values of each SNP with each predictor tested. Then we can assign the highest scoring predictor to the SNP. We then can save the table of candidates for downstream analysis.
Note that this assignation of a predictor to each SNP is only for visualization purposes as the process is not exactly fail-proof. Many SNPs have similar correlation values with different predictors, but only one will be assigned as correlated. When we cross the list of candidates with other uni-variate GEA analyses (BayPass), we can see in what univariate analysis each SNP was also recovered as a candidate.
```{R}
# add correlation of candidate SNPs with all environmental variables
foo <- matrix(nrow=(ncand), ncol=NCOL(env.predictors))
colnames(foo) <- colnames(env.predictors)
for (i in 1:length(cand$snp)) {
  nam <- gsub(":", ".", cand[i,2])
  snp.gen <- gt_data_tsv[,nam]
  foo[i,] <- apply(env.predictors, 2, function(x) cor(x, snp.gen))
}
cand <- cbind.data.frame(cand,foo)  

# assign SNPs to predictor based on highest correlation value - might be worth exploring
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,(4+(NCOL(env.predictors)))] <- names(which.max(abs(bar[4:(3+(NCOL(env.predictors)))]))) # gives the variable
  cand[i,(5+(NCOL(env.predictors)))] <- max(abs(bar[4:(3+(NCOL(env.predictors)))]))              # gives the correlation
}
colnames(cand)[(4+(NCOL(env.predictors)))] <- "predictor"
colnames(cand)[(5+(NCOL(env.predictors)))] <- "correlation"

# write table of candidate snps - 2 PCs
write.table(x = cand,
            file = "candidates_twopcs_table.tsv",
            quote=FALSE,  col.names = T, row.names = F, sep= "\t")

# write table of candidate snps - 5 PCs
write.table(x = cand,
            file = "candidates_fivepcs_table.tsv",
            quote=FALSE,  col.names = T, row.names = F, sep= "\t")

# write counts of candidates assigned to each variable - 2 PCs
write.table(x = table(cand$predictor),
            file = "number_of_candidates_twopcs_table.tsv",
            quote=FALSE,  col.names = F, row.names = F, sep= "\t")

# write counts of candidates assigned to each variable - 5 PCs
write.table(x = table(cand$predictor),
            file = "number_of_candidates_fivepcs_table.tsv",
            quote=FALSE,  col.names = F, row.names = F, sep= "\t")
```
We can transform our candidates table to a BED file for downstream analysis in bash
```{bash}
cut -f1,2 candidates_table.tsv | grep -vw "snp" | rev | cut -d'_' -f2- | rev | tr ':' '\t' | awk '{FS="\t"; OFS="\t"; print $2, $3-1, $3, $1;}' > rda_candidate_snps.bed
```
To plot RDA results we first need to prepare the data so we can discern between candidates and non-candidate SNPs and also assign a different color to each correlated predictor.
```{R}
# candidate snps
sel <- cand$snp
# predictors of each snp
env <- cand$predictor
# choose and add colors for each predictor
env[env=="PC1"] <- 'purple'
env[env=="PC2"] <- 'chartreuse3'
env[env=="PC3"] <- 'coral3'
env[env=="PC4"] <- 'gold'
env[env=="PC5"] <- 'dodgerblue4'

# pull all the SNP names
col.pred <- rownames(rda$CCA$v)
# replace snp name with its color code (only candidates here)
for (i in 1:length(sel)) {
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}
# replace non-candidates with grey
col.pred[grep("affold",col.pred)] <- '#f1eef6'

# create an empty (transparent) list for when drawing only candidates
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
# and same but with grey for candidates
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")

# list all colors 2 PCs
bg <- c('purple','chartreuse3')

# list all colors 5 PCs
bg <- c('purple','chartreuse3','coral3','gold','dodgerblue4')
```
Then we can plot the SNPs within RDA1 v RDA2 and RDA1 v RDA3 axes
```{R}
# axes 1 & 2 - 2 PCs
pdf(file = paste0("twopcs_rda1_rda2_snps.pdf"),
    width = 8,
    height = 8)

# axes 1 & 2 - 5 PCs
pdf(file = paste0("fivepcs_rda1_rda2_snps.pdf"),
    width = 8,
    height = 8)

plot(rda, type="n", 
     xlim=c(-(abs(load.rda[which.min(load.rda[,1]),1])*1.5),(abs(load.rda[which.max(load.rda[,1]),1])*1.5)), 
     ylim=c(-(abs(load.rda[which.min(load.rda[,2]),1])*1.5),(abs(load.rda[which.max(load.rda[,2]),1])*1.5)))
points(rda, display="species", pch=4, cex=0.7, col="gray32", bg=col.pred)
points(rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty)
text(rda, scaling=3, display="bp", col="#0868ac", cex=1)
# 2 PCs
legend("bottomright", legend=c("PC1", "PC2"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# 5 PCs
legend("bottomright", legend=c("PC1", "PC2", "PC3", "PC4", "PC5"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

dev.off()

# axes 1 & 3 - 5 PCs
pdf(file = paste0("fivepcs_rda1_rda3_snps.pdf"),
    width = 8,
    height = 8)


plot(rda, type="n", 
     xlim=c(-(abs(load.rda[which.min(load.rda[,1]),1])*1.5),(abs(load.rda[which.max(load.rda[,1]),1])*1.5)), 
     ylim=c(-(abs(load.rda[which.min(load.rda[,3]),1])*1.5),(abs(load.rda[which.max(load.rda[,3]),1])*1.5)), choices=c(1,3))
points(rda, display="species", pch=4, cex=0.7, col="gray32", bg=col.pred, choices=c(1,3))
points(rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, choices=c(1,3))
text(rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("PC1", "PC2", "PC3", "PC4", "PC5"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

dev.off()
```
We can transform our candidate table into a BED for future analyses
```{bash}
cut -f1,2 candidates_fivepcs_table.tsv | grep -vw "snp" | rev | cut -d'_' -f2- | rev |
 tr ':' '\t' | awk '{FS="\t"; OFS="\t"; print $2, $3-1, $3, $1;}' > rda_candidate_fivepcs_snps.bed
```
