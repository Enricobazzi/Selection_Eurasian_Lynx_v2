---
title: "GenWin"
author: "Enrico"
date: "12/27/2021"
output: html_document
editor_options:
  chunk_output_type: console
---
GenWin is a program that trys to detect inflection points in summary statistics calculated in a locus per locus manner (e.g. Fst). Using these inflection points it will divide the genome in windows which are not of arbitrary size, but based on properties of the summary statistic values. The program will output the number of SNPs inside each window and give a summary value based on the spline (Wstat). Wstat is a value that depends on the mean of the summary statistic, weighted against the mean and standard deviation of the dataset and the number of SNPs inside each window.

Wstat can be then used to calculate outlier windows with a quantile criteria.

I will use it to divide the Lynx lynx genome in windows based on the Bayes Factor (probability of a SNP being associated with a particular environmental predictor - See BayPass.md), calculated by BayPass, as a summary statistic.
```{R}
library(tidyverse)
library(GenWin)

# Define the variables to work in a loop
# WorldClim + Snow:
variables <- c("bio2", "bio5", "bio6", "bio8", "bio13", "jan_depth", "snow_days")

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
                          plotRaw=T, plotWindows=T,  method=3)

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
for var in bio2 bio5 bio6 bio8 bio13 jan_depth snow_days
 do
  echo "${var}"
  grep -v "WindowStart" ${var}_GenWin_windows_outliers.tsv |
   awk '{FS="\t"; OFS="\t"; print $6, $1, $2;}' > ${var}_GenWin_windows_outliers.bed
done
scp ~/Documents/Selection_Eurasian_Lynx_v2/3-Identify_Candidate_Loci/tables/*_GenWin_windows_outliers.bed \
ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/GenWin/
```
Extract the TOP SNP of each window (highest BF value), including only one SNP for the entire inversion (found between bases 17195000 and 18625000 on scaffold scaffold_17_arrow_ctg1).
```{R}
variables <- c("bio2", "bio5", "bio6", "bio8", "bio13", "jan_depth", "snow_days")

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

 # Results BED file:
 bedresults <- data.frame(chr=snps.table$scaffold, start=(snps.table$position -1),
                          stop=snps.table$position, BF=snps.table$BF.dB.)
 
 outlier_windows <- read.table(paste0("3-Identify_Candidate_Loci/tables/",var,"_GenWin_windows_outliers.tsv"),
                               h=T, as.is = T)
 
 topsnps <- data.frame()
 for (i in 1:nrow(outlier_windows)){
   window <- outlier_windows[i,]
   scaffold <- window$scaffold
   Wstart <- window$WindowStart
   Wstop <- window$WindowStop
   if (scaffold=="scaffold_17_arrow_ctg1" & Wstart >= 17195000 & Wstop <= 18625000) {
     window_snps <- subset(bedresults, chr == scaffold & start >= 17195000 & stop <= 18625000)
   } else {
   window_snps <- subset(bedresults, chr == scaffold & start >= Wstart & stop <= Wstop)
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
scp ~/Documents/Selection_Eurasian_Lynx_v2/3-Identify_Candidate_Loci/tables/*_topsnps.range \
ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/GenWin/
```

