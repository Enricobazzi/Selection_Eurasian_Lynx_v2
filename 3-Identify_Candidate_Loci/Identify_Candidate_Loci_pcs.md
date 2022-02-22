---
title: "Identify_Candidate_Loci_pcs"
author: "Enrico"
date: "2/18/2022"
output: html_document
editor_options:
  chunk_output_type: console
---

The strategy adopted to identify loci under selection is a Genome Environment Association (GEA) analysis using the first 5 axes of the PCA of environmental variables as independent predictors.

## Redundancy Analysis (RDA)

Of the different methods reviewed to perform GEA, it appears that redundancy analysis (RDA) is one of the most reliable (Forester et al. 2018). RDA is a multivariate method for exploring association of different environmental predictors with genetic data. By exploring which SNPs more strongly correlate with the different constrained axes of RDA ordination, we can identify candidates for selection.

I followed a few tutorials made to guide through the different steps of GEA and RDA analysis. Links of the tutorials are:

https://popgen.nescent.org/2018-03-27_RDA_GEA.html
https://bookdown.org/hhwagner1/LandGenCourse_book/WE-11.html

See RDA_pcs.md for details on how the analysis was run and its outputs.

## Univariate GEA with BayPass

In order to reduce false positive rates and add a control for population structure (RDA does not control for it), I have run univariate analysis using the software BayPass v2.1. See BayPass_pcs.md for details on how the analysis was run and its outputs.

BayPass results, in the form of SNPs Bayes Factors, have been used to generate a set of candidate genomic windows, using the software GenWin. See GenWin_pcs.md for details on how the analysis was run and its outputs.

## Intersecting results to obtain candidate windows and SNPs

By intersecting our RDA results with the genomic windows identified by BayPass, we can point which genomics windows and which SNPs are under selection and run downstream analysis to answer different questions. What genes are being selected and what function do they have for the organism? What effect does the environment have on the genetic structure of the populations? Can we identify an environmental gradient that is directly related to the genetic structure of the populations?

```{bash}
# On genomics-a server
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/Intersect

for var in PC1 PC2 PC3 PC4 PC5
 do
  bedtools intersect -a ../RDA/rda_candidate_fivepcs_snps.bed -b ../GenWin/${var}_GenWin_windows_outliers.bed \
   > ${var}_intersect_candidate_fivepcs_snps.bed
  bedtools intersect -wa -a ../GenWin/${var}_GenWin_windows_outliers.bed -b ../RDA/rda_candidate_fivepcs_snps.bed | uniq \
   > ${var}_intersect_candidate_fivepcs_windows.bed
  nsnps=($(wc -l <${var}_intersect_candidate_fivepcs_snps.bed))
  nwindows=($(wc -l <${var}_intersect_candidate_fivepcs_windows.bed))
  winlength=($(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${var}_intersect_candidate_fivepcs_windows.bed))
  echo "${var} has a total of ${nsnps} candidate SNPs, a total of ${nwindows} candidate genomic windows \
  spanning a total of ${winlength} bps"
done
cat *_intersect_candidate_fivepcs_snps.bed | sort -k 1,1 -k2,2n | uniq \
 > total_intersect_candidate_fivepcs_snps.bed
cat *_intersect_candidate_fivepcs_windows.bed | sort -k 1,1 -k2,2n | uniq | bedtools merge -i - \
 > total_intersect_candidate_fivepcs_windows.bed
```
PC1 has a total of 1813 candidate SNPs, a total of 320 candidate genomic windows spanning a total of 7480000 bps
PC2 has a total of 1547 candidate SNPs, a total of 334 candidate genomic windows spanning a total of 7500000 bps
PC3 has a total of 2701 candidate SNPs, a total of 484 candidate genomic windows spanning a total of 11320000 bps
PC4 has a total of 872 candidate SNPs, a total of 268 candidate genomic windows spanning a total of 6160000 bps
PC5 has a total of 5217 candidate SNPs, a total of 715 candidate genomic windows spanning a total of 15040000 bps

total candidate SNPs : 9871
total candidate genomic windows : 1427
total length of candidate genomic windows : 42.27 Mbp

## Get the TOPSNP of each candidate window

Top SNPs for candidate windows of each PC were extracted (see GenWin_pcs.md). To put together all of them into a single file:
```{bash}
# on genomics-a:
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/

# get topsnps from all uncorrelated predictors
rm Intersect/fivepcs_topsnps.range
for var in PC1 PC2 PC3 PC4 PC5
 do
  echo "${var}"
  cat GenWin/${var}_topsnps.range >> Intersect/fivepcs_topsnps.range
done
```
Problem = there are more than one "topsnp" for each candidate window (either consecutive windows with different "topsnps" were joined, or the window is the same but the "topsnp" is different for different PCs)
Solution = get topsnps from each candidate window with only one topnsp and then get one random topsnp from each window with more than one

First we have to join all the windows from the inversion (see GenWin_pcs.md to see how they were defined):
between 17635000 and 18355000 of scaffold "scaffold_17_arrow_ctg1"

I manually modified "total_intersect_candidate_fivepcs_windows.bed" to have a single window for the inversion using nano - resulting in a total of 1422 unique windows down from 1427 (6 windows in the inversion)

Then I can run:
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/Intersect

# get list of candidate windows with ONE topsnp
bedtools intersect -wo -a total_intersect_candidate_fivepcs_windows.bed \
 -b <(cat fivepcs_topsnps.range | sort -k 1,1 -k2,2n) |
 cut -f1,2,3 | uniq -u > topsnps_fivepcs_nodupwindows.bed
 
# get topsnp from each window with ONE topsnp
bedtools intersect -a <(cat fivepcs_topsnps.range | sort -k 1,1 -k2,2n) \
 -b topsnps_fivepcs_nodupwindows.bed > topsnps_fivepcs_nodupwindows_onesnp.bed

# get list of candidate windows with MORE than one topsnp
bedtools intersect -wo -a total_intersect_candidate_fivepcs_windows.bed \
 -b <(cat fivepcs_topsnps.range | sort -k 1,1 -k2,2n) |
 cut -f1,2,3 | uniq -d > topsnps_fivepcs_dupwindows.bed

# get only one topsnp for each of the windows with duplicates
while read p; do
  bedtools intersect -a <(cat fivepcs_topsnps.range | sort -k 1,1 -k2,2n) -b <(echo "$p") | shuf -n 1
done < topsnps_fivepcs_dupwindows.bed > topsnps_fivepcs_dupwindows_onerand.bed

# get full list of topsnps (only one per window)
cat topsnps_fivepcs_nodupwindows_onesnp.bed topsnps_fivepcs_dupwindows_onerand.bed | sort -k 1,1 -k2,2n \
 > total_intersect_candidate_fivepcs_windows_topsnps.bed
```
Total of 1422 topsnps found from 1422 unique candidate windows

Download results to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/Intersect/total_intersect_candidate_fivepcs_windows.bed ~/Documents/Selection_Eurasian_Lynx_v2/3-Identify_Candidate_Loci/tables/
```

## Plotting results

```{R}
library(tidyverse)
library(venn)

variables <- c("PC1", "PC2", "PC3", "PC4", "PC5")

## Manhattan plotting of Bayes Factors of SNPs

for (k in 1:length(variables)){

 # Variable:
 var=variables[k]

 # Empty SNPs dataframe to fill
 snps.table <- data.frame()
 
 # Fill SNPs dataframe with all BayPass results for the variable
 for (n in 1:50){
  # upload the dataset table:
  var.snp=read.table(paste0("3-Identify_Candidate_Loci/tables/BayPass_OutPut/AUX_",
                            var,"_",n,"_summary_betai.out"),h=T)
  
  # calculate the SNP database number sequence (1 every 50) 
  # NOTE: 2100553 is the total number of SNPs in all datasets
  var.snp <- data.frame(var.snp, SNPnum = seq(n, 2100553, by = 50))
  
  # add the dataset rows to the total snps table
  snps.table <- rbind(snps.table, var.snp)
 }
 # order the table based on the SNP number
 snps.table <- snps.table %>% arrange(SNPnum)
 
 # Add SNP ID information from SNPIDs table
 SNPIDs <- read_tsv("3-Identify_Candidate_Loci/tables/BayPass_OutPut/finalset.maf5pc.SNPIDs",
                    col_names = F)[,2-3] %>%
    rename("scaffold" =  X2, "position" = X3)
  
 # Add SNP IDs to the total snps table
 snps.table <- data.frame(snps.table,SNPIDs)

 # plot
 p <- ggplot() +
   geom_point(data=snps.table, aes(x=SNPnum, y=BF.dB.), fill="grey32", shape=21, size=1) +
   theme_minimal()
 
 ggsave(p, filename = paste0("3-Identify_Candidate_Loci/plots/BFs_", var, "_manplot.pdf"),  width = 14,
     height = 4)
 
}


## Manhattan plotting outlier and candidate windows ##

for (k in 1:length(variables)){

 # Variable:
 var=variables[k]

 # Windows
 windows <- read.table(paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows.tsv"),
                               h=T, as.is = T)
 
 # Outliers
 outlier_windows <- read.table(paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows_outliers.tsv"),
                               h=T, as.is = T)
 
 # Outliers intersected
 intersect_table <- read.table(paste0("3-Identify_Candidate_Loci/tables/total_intersect_candidate_windows.bed"),
                               h=F, as.is = T)
 
 intersect_windows <- data.frame()
 for (n in 1:NROW(intersect_table)){
   row <- intersect_table[n,]
   SUB <- subset(outlier_windows, scaffold == row$V1 & WindowStart >= row$V2 & WindowStop <= row$V3)
   intersect_windows <- rbind(intersect_windows, SUB)
 }
 
 # plot
 p <- ggplot() +
   geom_point(data=windows, aes(x=WindowNumber, y=Wstat), color="grey32", shape=4, size=1.5) +
   geom_point(data=outlier_windows, aes(x=WindowNumber, y=Wstat), color="gold2",shape=4, size=1.5) +
   geom_point(data=intersect_windows, aes(x=WindowNumber, y=Wstat), color="brown", fill="gold2", shape=21, size=3) +
   theme_minimal()
 
 ggsave(p, filename = paste0("3-Identify_Candidate_Loci/plots/outliers_intersect_", var, "_manplot.pdf"),  width = 14,
     height = 4)
 
}

## Plot intersections with Venn Diagrams ##

variables <- c("PC1", "PC2", "PC3", "PC4", "PC5")

v.list <- list()

for (k in 1:length(variables)){
 
 # Variable:
 var=variables[k]

 # Outliers
 outlier_windows <- read.table(paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows_outliers.tsv"),
                               h=T, as.is = T)
 
 # Outliers intersected
 intersect_table <- read.table(paste0("3-Identify_Candidate_Loci/tables/total_intersect_candidate_windows.bed"),
                               h=F, as.is = T)
 
 intersect_table$V4 <- 1:NROW(intersect_table)
 
 intersect_windows <- data.frame()
 for (n in 1:NROW(outlier_windows)){
  row <- outlier_windows[n,]
  SUB <- subset(intersect_table, V1 == row$scaffold & V2 <= row$WindowStart & V3 >= row$WindowStop)
  intersect_windows <- rbind(intersect_windows, SUB)
 }

 # Create vector of intersect window numbers
 vec <- c(as.character(intersect_windows$V4))
 
 # add vector to list
 v.list <- c(v.list, list(vec))
}
length(unique(unlist(v.list)))

pdf(file = paste0("3-Identify_Candidate_Loci/plots/intersect_fivepcs_venn_diagram.pdf"),
    width = 8,
    height = 8)
venn(v.list, snames = variables, 
     #zcolor = "lightseagreen, chartreuse3, coral, firebrick3, purple, gold, dodgerblue4",
     zcolor = 'style',
     ilcs = 0.9,
     box = F)
dev.off()
```

## Plotting INVERSION

```{R}
library(tidyverse)

## For BayPass results:

variables <- c("PC1", "PC2", "PC3", "PC4", "PC5")

for (k in 1:length(variables)){
 # Variable:
 var=variables[k]

 # Windows
 snps.table <- data.frame()
 
 # Fill SNPs dataframe with all BayPass results for the variable
 for (n in 1:50){
  # upload the dataset table:
  print(n)
  var.snp=read.table(paste0("3-Identify_Candidate_Loci/tables/BayPass_OutPut_pcs/AUX_",
                            var,"_",n,"_summary_betai.out"),h=T)
  
  # calculate the SNP database number sequence (1 every 50) 
  # NOTE: 2100553 is the total number of SNPs in all datasets
  var.snp <- data.frame(var.snp, SNPnum = seq(n, 2100553, by = 50))
  
  # add the dataset rows to the total snps table
  snps.table <- rbind(snps.table, var.snp)
 }
 # order the table based on the SNP number
 snps.table <- snps.table %>% arrange(SNPnum)
 
 # Add SNP ID information from SNPIDs table
 SNPIDs <- read_tsv("3-Identify_Candidate_Loci/tables/BayPass_OutPut/finalset.maf5pc.SNPIDs",
                    col_names = F)[,2-3] %>%
    rename("scaffold" =  X2, "position" = X3)
  
 # Add SNP IDs to the total snps table
 snps.table <- data.frame(snps.table,SNPIDs)
 
 inversion_snps <- subset(snps.table, scaffold == "scaffold_17_arrow_ctg1" & position >= 14005000 & position <= 22000000)

 p <- ggplot() +
   geom_point(data=inversion_snps, aes(x=position, y=BF.dB.), fill="grey32", shape=21, size=1.5) +
   theme_minimal()
 ggsave(p, filename = paste0("3-Identify_Candidate_Loci/plots/inversion_", var, "_manplot.pdf"),  width = 14,
     height = 4)
}

## For RDA results

rda.loadings <- read_tsv("3-Identify_Candidate_Loci/tables/rda_loadings_sigaxes.tsv",
                    col_names = T)

rda.loadings.final <- data.frame()

for (n in 1:length(rda.loadings$SNP)){
  print(n)
  snp.row <- rda.loadings[n,]
  snp <- rda.loadings$SNP[n]
  SCAF <- strsplit(snp, split = ":")[[1]][1]
  POS <- (strsplit(strsplit(snp, split = ":")[[1]][2], split = "_"))[[1]][1]
  row <- data.frame(scaffold=SCAF, position=POS, RDA1=snp.row$RDA1, RDA2=snp.row$RDA2, RDA3=snp.row$RDA3)
  rda.loadings.final <- rbind(rda.loadings.final, row)
}

```
