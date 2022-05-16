---
title: "GenWin_pcs"
author: "Enrico"
date: "2/18/2022"
output: html_document
editor_options:
  chunk_output_type: console
---
## Introduction to GenWin
GenWin is a program that trys to detect inflection points in summary statistics calculated in a locus per locus manner (e.g. Fst). Using these inflection points it will divide the genome in windows which are not of arbitrary size, but based on properties of the summary statistic values. The program will output the number of SNPs inside each window and give a summary value based on the spline (Wstat). Wstat is a value that depends on the mean of the summary statistic, weighted against the mean and standard deviation of the dataset and the number of SNPs inside each window.

Wstat can be then used to calculate outlier windows with a quantile criteria.

I will use it to divide the Lynx lynx genome in windows based on the Bayes Factor (probability of a SNP being associated with a particular environmental predictor - See BayPass.md), calculated by BayPass, as a summary statistic.

## Running GenWin

```{R}
library(tidyverse)
library(GenWin)

# Define the variables to work in a loop
# WorldClim + Snow:
variables <- c("PC1", "PC2", "PC3", "PC4", "PC5")

for (k in 1:length(variables)){

 # Variable:
 var=variables[k]
 print(var)
 
 # Empty SNPs dataframe to fill
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
 
 # Empty Spline output table to fill with GenWin results
 all_spline <- data.frame()
 # Empty Spline output table to fill with the OUTLIERS of GenWin results
 all_spline_outliers <- data.frame()
 
 # Loop through all scaffolds as GenWin works on one scaffold at the time
 for (n in 1:length(unique(snps.table$scaffold))) {
 
  # Scaffold:
  chr=unique(snps.table$scaffold)[n]
  print(chr)
  
  # Get SNPs table for the scaffold only
  data_chr <- snps.table %>% filter(scaffold==chr)
 
  # Analyze data with GenWin - Produces a table and a plot
  spline <- splineAnalyze(data_chr$BF.dB., data_chr$position,
                          smoothness=10000, 
                          mean=mean(snps.table$BF.dB.),
                          s2=var(snps.table$BF.dB.),
                          method=3)

  # Get the table with the results
  spline.data <- as.data.frame(spline$windowData)
  # Add the scaffold name to the results table
  spline.data$scaffold <- rep(chr,nrow(spline.data))
  # Add the variable
  spline.data$var <- rep(var,nrow(spline.data))
  # Add the Window length
  spline.data$WindowLength <- spline.data$WindowStop-spline.data$WindowStart

  # Alternate Odd and Even - for later plot of all scaffolds
  if((n %% 2) == 0) {
   spline.data$color <- rep("Even",nrow(spline.data))
  } else {
   spline.data$color <- rep("Odd",nrow(spline.data))
  } 

  # Fill the GenWin results table
  all_spline <- data.frame(rbind(all_spline, data.frame(spline.data)))
 }
 
 # Add the Window Number to spline dataset
 all_spline$WindowNumber <- 1:nrow(all_spline)
 
 # Calculate threshold for all genome
 all_spline.thresh=as.numeric(quantile(all_spline$Wstat,probs=0.99,na.rm=T))
 # Get whole-genome outliers table
 all_spline_total_outliers <- subset(all_spline, Wstat > all_spline.thresh)
 
 # Save the tables generated
 # All windows:
 write.table(x = all_spline,file = paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")
 # Whole-Genome threshold outlier windows:
 write.table(x = all_spline_total_outliers,
             file = paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows_outliers.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")
}
```
Convert to BED and upload to Genomics server
```{bash}
cd /Users/enrico/Documents/Selection_Eurasian_Lynx_v2/3-Identify_Candidate_Loci/tables

for var in PC1 PC2 PC3 PC4 PC5
 do
  echo "${var}"
  grep -v "WindowStart" ${var}_GenWin_windows_outliers.tsv |
   awk '{FS="\t"; OFS="\t"; print $6, $1, $2;}' > ${var}_GenWin_windows_outliers.bed
done
scp ~/Documents/Selection_Eurasian_Lynx_v2/3-Identify_Candidate_Loci/tables/*_GenWin_windows_outliers.bed \
ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/GenWin/
```

## Define inversion Windows

```{R}
## For GenWin windows of BayPass results:
# Explore the area around the inversion (from 14005000 to 220000000 of scaffold_17_arrow_ctg1)
k=5
for (k in 1:length(variables)){
 # Variable:
 var=variables[k]
 
 # Explore manually outlier windows for consecutive outliers
 windows_outliers <- read.table(paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows_outliers.tsv"),
                               h=T, as.is = T)
 inversion_windows_outliers <- subset(windows_outliers,
                             scaffold == "scaffold_17_arrow_ctg1" & WindowStart >= 14005000 & WindowStop <= 22000000)
 
 # Plot all windows in area for visual exploration
 windows <- read.table(paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows_outliers.tsv"),
                               h=T, as.is = T)
 inversion_windows <- subset(windows,
                             scaffold == "scaffold_17_arrow_ctg1" & WindowStart >= 14005000 & WindowStop <= 22000000)
  p <- ggplot() +
   geom_point(data=inversion_windows, aes(x=WindowStart, y=Wstat), fill="grey32", shape=21, size=1.5) +
   theme_minimal()

 ggsave(p, filename = paste0("3-Identify_Candidate_Loci/plots/inversion_", var, "_windows_manplot.pdf"),  width = 14,
     height = 4)
}
# PC1 17635000 to 18355000 + 18415000 to 18445000 + 18615000 to 18625000
# PC2 17635000 to 18355000 + 18615000 to 18625000
# PC3 18405000 to 18435000 - outside of window?
# PC4 17635000 to 18365000 - more non-outlier windows in the middle of this range + 18615000 to 18625000
# PC5 17465000 to 17535000 as outliers (before rest) - but whole area higher but not outlier
```

I define as the inversion the range between 17635000 and 18355000 of scaffold "scaffold_17_arrow_ctg1"

To verify if it's really an inversion, I made a PCA of the SNPs present in this range. See PCA_inversion.R

## Extract TOP SNP of each Window

Extract the TOP SNP of each window (highest BF value), including only one SNP for the entire inversion (found between bases 17635000 and 18355000 on scaffold scaffold_17_arrow_ctg1) and excluding any SNPs with missing data.

To get list of SNPs with missing data
```{bash}
# on genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF
# remove samples from ll_wholegenome_LyCa_ref.sorted.filter7.vcf
samplesARRAY=($(grep -m1 "#CHROM" ll_wholegenome_LyCa_ref.sorted.filter7.vcf | tr '\t' '\n' | grep "_ll_" | grep -vE "h_ll_ba_0214|h_ll_ba_0215|c_ll_ba_0216"))

/opt/gatk-4.1.0.0/gatk SelectVariants \
-R /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa \
-V ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
$(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
-O all_samples.vcf


# get SNPs from the finalset.maf5pc from allsamples vcf
bedtools intersect -a all_samples.vcf \
 -b ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf \
 > all_samples_finalsetvariants.vcf

# get column defining where is missing data  
grep -v "#" all_samples_finalsetvariants.vcf | cut -d';' -f3 | cut -d'=' -f2 | awk '{ if ($0 == "210") print "all"; else print "missing"; }' > missing_data_col.txt
```
In R get top snps:
```{R}
variables <- c("PC1", "PC2", "PC3", "PC4", "PC5")

for (k in 1:length(variables)){

 # Variable:
 var=variables[k]
 # Empty SNPs dataframe to fill
 snps.table <- data.frame()
 
 # Fill SNPs dataframe with all BayPass results for the variable
 for (n in 1:50){
  # upload the dataset table:
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
    rename(X2 = "scaffold", X3 = "position")
  
 # Add column for missing data
 miss_data_col <- read_tsv("3-Identify_Candidate_Loci/tables/missing_data_col.txt",
                    col_names = F)
 # Add SNP IDs to the total snps table
 snps.table <- data.frame(snps.table,SNPIDs,miss_data_col)

 # Results BED file:
 bedresults <- data.frame(chr=snps.table$scaffold, start=(snps.table$position -1),
                          stop=snps.table$position, BF=snps.table$BF.dB., miss=snps.table$X1)
 
 outlier_windows <- read.table(paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows_outliers.tsv"),
                               h=T, as.is = T)
 
 topsnps <- data.frame()
 for (i in 1:nrow(outlier_windows)){
   window <- outlier_windows[i,]
   scaffold <- window$scaffold
   Wstart <- window$WindowStart
   Wstop <- window$WindowStop
   if (scaffold=="scaffold_17_arrow_ctg1" & Wstart >= 17635000 & Wstop <= 18355000) {
     window_snps <- subset(bedresults, chr == scaffold & start >= 17635000 & stop <= 18355000)
   } else {
   window_snps <- subset(bedresults, chr == scaffold & start >= Wstart & stop <= Wstop & miss == "all")
   }
   topsnp <- window_snps[order(window_snps$BF,decreasing = TRUE),]
   topsnp <- topsnp[1,]
   topsnp$var <- var
   topsnp$BF <- NULL
   topsnps <- rbind(topsnps, topsnp)
 }
 topsnps <- distinct(topsnps)
 write.table(x = topsnps,
             file = paste0("3-Identify_Candidate_Loci/tables/",var,"_topsnps.range"),
             quote=FALSE,  col.names = F, row.names = FALSE, sep= "\t")
}
```
Upload to Genomics server
```{bash}
scp ~/Documents/Selection_Eurasian_Lynx_v2/3-Identify_Candidate_Loci/tables/PC*_topsnps.range \
ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/GenWin/
```
