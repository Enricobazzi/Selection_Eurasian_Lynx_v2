---
title: "Downstream_Analyses"
author: "Enrico"
date: "7/12/2021"
output: html_document
editor_options:
  chunk_output_type: console
---

Downstream analyses of the candidate loci include different approaches to explore mainly 3 different aspects:
  
  1. Effects of selection on the genetic structure of populations
  2. Spatial gradients of environmental effects on genetics
  3. Functional annotation of the candidate loci

## 1. Effects of selection on the genetic structure of populations

To explore this aspect we will run RDA again, but this time we will compare structure recovered from selected vs unselected loci.

## RDA on candidate SNPs only

```{bash}
# on genomics-a:
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/

# get topsnps from all uncorrelated predictors
for var in bio2 bio5 bio6 bio8 bio13 jan_depth snow_days
 do
  echo "${var}"
  cat GenWin/${var}_topsnps.range >> Intersect/uncorvars_topsnps.range
done
```
Problem = more than one topsnp for each candidate window (window is same topsnp is different because from different variable or originally window was split in 2)
Solution = get topsnps from each candidate window with only one topnsp and then get one random topsnp from each window with more than one
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/Intersect

# get list of candidate windows with ONE topsnp
bedtools intersect -wo -a total_intersect_candidate_windows.bed \
 -b <(cat uncorvars_topsnps.range | sort -k 1,1 -k2,2n) |
 cut -f1,2,3 | uniq -u > topsnps_nodupwindows.bed
 
# get topsnp from each window with ONE topsnp
bedtools intersect -a <(cat uncorvars_topsnps.range | sort -k 1,1 -k2,2n) \
 -b topsnps_nodupwindows.bed > topsnps_nodupwindows_onesnp.bed

# get list of candidate windows with MORE than one topsnp
bedtools intersect -wo -a total_intersect_candidate_windows.bed \
 -b <(cat uncorvars_topsnps.range | sort -k 1,1 -k2,2n) |
 cut -f1,2,3 | uniq -d > topsnps_dupwindows.bed

# get only one topsnp for each of the windows with duplicates
while read p; do
  bedtools intersect -a <(cat uncorvars_topsnps.range | sort -k 1,1 -k2,2n) -b <(echo "$p") | shuf -n 1
done < topsnps_dupwindows.bed > topsnps_dupwindows_onerand.bed

# get full list of topsnps (only one per window)
cat topsnps_nodupwindows_onesnp.bed topsnps_dupwindows_onerand.bed | sort -k 1,1 -k2,2n \
 > total_intersect_candidate_windows_topsnps.bed
```
Total of 211 topsnps found from 211 unique candidate windows - although initially 220 candidate windows, 9 of them are from the inversion so only 1 SNP will be extracted from them

Extract VCF and PLINK of the topsnps
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/

bedtools intersect -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf -b Intersect/total_intersect_candidate_windows_topsnps.bed \
 > VCF/ll_allsamples_total_intersect_candidate_windows_topsnps.vcf

cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/
# I manually created a file listing samples to remove for plink with following format:
# c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

plink_1.9 --vcf ll_allsamples_total_intersect_candidate_windows_topsnps.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--remove samplestoremove.txt --geno 0 \
--recode A --out total_intersect_candidate_windows_topsnps
```
Copy to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/total_intersect_candidate_windows_topsnps.raw Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
Prepare R for RDA on laptop
```{R}
# load libraries
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(adegenet)
library(PopGenReport)
library(gdistance)

# Load environmental data
env.predictors <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")

# Load GenoType Data - choose set of vars file
gt_data <- read.PLINK("4-Downstream_Analyses/tables/total_intersect_candidate_windows_topsnps.raw")
gt_data_tsv <- data.frame(as.matrix(gt_data))
```
Run RDA
```{R}
rda <- rda(gt_data ~ ., data=env.predictors, scale=T)
signif.axis <- anova.cca(rda, by="axis", parallel=getOption("mc.cores"))
```
Prepare plot
```{R}
# Add populations and colors
loc <- rep(NA, NROW(gt_data_tsv))
loc[grep("ba", rownames(gt_data_tsv))] <- "Balkans"
loc[grep("ca", rownames(gt_data_tsv))] <- "Caucasus"
loc[grep("cr", rownames(gt_data_tsv))] <- "Carpathians"
loc[grep("ka", rownames(gt_data_tsv))] <- "Mongolia"
loc[grep("ki", rownames(gt_data_tsv))] <- "Kirov"
loc[grep("la", rownames(gt_data_tsv))] <- "Latvia"
loc[grep("no", rownames(gt_data_tsv))] <- "Norway"
loc[grep("po", rownames(gt_data_tsv))] <- "NE-Poland"
loc[grep("og", rownames(gt_data_tsv))] <- "Mongolia"
loc[grep("to", rownames(gt_data_tsv))] <- "Mongolia"
loc[grep("tu", rownames(gt_data_tsv))] <- "Tuva"
loc[grep("ur", rownames(gt_data_tsv))] <- "Urals"
loc[grep("vl", rownames(gt_data_tsv))] <- "Vladivostok"
loc[grep("ya", rownames(gt_data_tsv))] <- "Yakutia"

num <- rep(NA, NROW(gt_data_tsv))
num[grep("ba", rownames(gt_data_tsv))] <- 1
num[grep("ca", rownames(gt_data_tsv))] <- 3
num[grep("cr", rownames(gt_data_tsv))] <- 2
num[grep("ka", rownames(gt_data_tsv))] <- 6
num[grep("ki", rownames(gt_data_tsv))] <- 4
num[grep("la", rownames(gt_data_tsv))] <- 5
num[grep("no", rownames(gt_data_tsv))] <- 8
num[grep("po", rownames(gt_data_tsv))] <- 7
num[grep("og", rownames(gt_data_tsv))] <- 6
num[grep("to", rownames(gt_data_tsv))] <- 6
num[grep("tu", rownames(gt_data_tsv))] <- 9
num[grep("ur", rownames(gt_data_tsv))] <- 10
num[grep("vl", rownames(gt_data_tsv))] <- 11
num[grep("ya", rownames(gt_data_tsv))] <- 12

ola <- data.frame(sample = rownames(gt_data_tsv), pop = loc, n = num)
eco <- levels(ola$pop)
bg <- cols <- c("#A035AF",
                brewer.pal(12,"Paired")[9],
                "#B8860b",
                viridis_pal()(5)[1],
                brewer.pal(12,"Paired")[3],
                brewer.pal(12,"Paired")[7],
                viridis_pal()(5)[3],
                viridis_pal()(5)[2],
                brewer.pal(12,"Paired")[8],
                "#0F4909",
                brewer.pal(12,"Paired")[5],
                brewer.pal(12,"Paired")[6])
```
plot
```{R}
# RDA1 v RDA2
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num]) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1)     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

# RDA1 v RDA3
plot(rda, type="n", scaling=3, choices=c(1,3))
points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3, choices=c(1,3))           # the SNPs
points(rda, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num], choices=c(1,3)) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

# RDA2 v RDA3
plot(rda, type="n", scaling=3, choices=c(2,3))
points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3, choices=c(2,3))           # the SNPs
points(rda, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num], choices=c(2,3)) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

# RDA3 v RDA4
plot(rda, type="n", scaling=3, choices=c(3,4))
points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3, choices=c(3,4))           # the SNPs
points(rda, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num], choices=c(3,4)) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(3,4))     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

# RDA5 v RDA6
plot(rda, type="n", scaling=3, choices=c(5,6))
points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3, choices=c(5,6))           # the SNPs
points(rda, display="sites", pch=21, cex=1.2, col="gray32", scaling=3, bg=bg[num], choices=c(5,6)) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(5,6))     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

```
