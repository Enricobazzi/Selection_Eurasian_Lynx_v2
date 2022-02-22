---
title: "BayPass_pcs"
author: "Enrico"
date: "2/10/2022"
output: html_document
editor_options:
  chunk_output_type: console
---
## Introduction to BayPass

In this Markdown I will describe and perform all the steps to perform a genome-wide scan for signatures of selection in Eurasian Lynx populations, using the software BayPass 2.1

The analyses will be performed on the genomics-b server of EBD, unless specified otherwise.

This software uses Allele Frequency data and Bayesian Hierarchical Models to generate a distribution of differentiation coefficients of SNPs across the genome. First a scaled population covariance matrix based on population allele frequencies is generated, to evaluate the relationships between populations. Based on this covariance matrix, and the supposed ancestral allele frequency (inferred from weighted mean of the reference allele frequency), a CORE model is generated.

SNPs that present values of differentiation (XtX, a SNP-specific Fst explicitly corrected for the scaled covariance of population allele frequencies) that exceed the amount expected under the core model can be identified as candidate selection loci.

Furthermore, this software allows the evaluation of the association of particular SNPs to some environmental Covariate. With one measurement for population for each covariate, the software will evaluate the data under to additional models:

(1) The Standard Covariate model (STD): which adds an association covariable (given by the correlation coefficient between the covariate measurement and XtX) to the CORE model.
(2) The Auxiliary Variable Covariate model (AUX): which further builds on the STD model by attaching a binary variable (0 or 1 -> association or no association) to each locus' regression coefficient. The posterior mean of this variable will indicate the posterior probability of the association of that variable with a particular SNP.

In order to run the software we will need:

– Allele Count data for all the considered populations: in the form of a space delimited file, with one row for each SNP and two columns for each population, with the allele counts for the reference and alternative alleles.

– Environmental variable files: in the form of different files with one column per population and their respective measurement for that particular covariate.

The Allele Count file will have to be generated from the VCF file containing information on all the populations of interest. Before generating the Allele Count file I will also proceed to filter out the samples I won't include in the final analysis, the regions that are not autosomic and variants with MAF < 5%, to avoid introducing noise from low frequency variants.

The Environmental variable files have already been generated and can be found at /home/ebazzicalupo/BayPass/Covariate_Data in the genomics-b server of EBD.

## Filtering for wanted samples, autosomic scaffolds and MAF > 5%

I will filter the final VCF from my variant filtering pipeline (2.Variant_filtering.md) to remove variants with a MAF of less than 5%. First I will subset my VCF to include only the samples I want and autosomic scaffolds only using GATK SelectVariants, then I will filter MAF using bcftools view:
```{bash}
cd /home/ebazzicalupo/BayPass/VCF

popARRAY=($(grep -m1 "#CHROM" ll_wholegenome_LyCa_ref.sorted.filter7.vcf | tr '\t' '\n' | grep "c_ll" | cut -d"_" -f3 | sort -u | grep -vE "ba|og|no|po|cr"))
samplesARRAY=($(for pop in ${popARRAY[@]}; do grep -m1 "#CHROM" ll_wholegenome_LyCa_ref.sorted.filter7.vcf | tr '\t' '\n' | grep "c_ll_${pop}" | grep -vE "c_ll_vl_0137|c_ll_tu_0154"; done))

/opt/gatk-4.1.0.0/gatk IndexFeatureFile -F ll_wholegenome_LyCa_ref.sorted.filter7.vcf

/opt/gatk-4.1.0.0/gatk SelectVariants \
-R /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa \
-V ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
$(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
-L /GRUPOS/grupolince/reference_genomes/lynx_canadensis/autosomic_scaffolds.bed \
-O ll_wholegenome_LyCa_ref.sorted.filter7.finalset.vcf


bcftools view -i 'MAF>0.05' ll_wholegenome_LyCa_ref.sorted.filter7.finalset.vcf \
> ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf
```

## Generating Allele Count data file

In order to generate the Allele Count file I will first have to divide my VCF with all the individuals, into different VCF for each population. Then the number of reference and alternative alleles can be computed and joined into a final file to be used as input for BayPass.

The script to perform this step can be found at Baypass_executables/Allele_Count_generation.sh and was launched as such

```{bash}
screen -S allelecounts
script allelecounts_baypass.log

./Allele_Count_generation.sh \
/home/ebazzicalupo/BayPass/VCF/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf \
/home/ebazzicalupo/BayPass/AlleleCounts
```

To avoid possible bias due to consecutive SNPs being non-independant observations because of Linkage Disequilibrium, I will divide my dataseta into 50 different dataset, made of one SNP every 50 (dataset1:snp1,snp51,snp101...;dataste2:snp2,snp52,snp102...;...).

To obtain these datasets I used the following script:

```{bash}
cd /home/ebazzicalupo/BayPass/AlleleCounts
# For the awk script to work I need to iterate from 0 to 49
for n in {0..49}
 do
  echo ${n}
  awk -v number="$n" 'NR % 50 == 0+number' all.allelecounts > all.allelecounts.${n}
done

# But the 0 iteration actually starts from the 50th line of the file
# so I will change the name of the file:
mv all.allelecounts.0 all.allelecounts.50
```

## CORE Model

Now I will run BayPass with the generated allele count data using the CORE model with no Co-Variate data.

```{bash}
# The total should be 1..50 - I will run it differently based on how many cores
# are available and the progress I made
for n in {1..50}
 do
  screen -dmS core_${n}  sh -c "/home/ebazzicalupo/BayPass/baypass_core.sh ${n}; exec /bin/bash"
done
```
The results under the CORE model can be analyzed in R on my laptop. The results directory (OutPut) was copied in the project directory for analysis with R.
```{R}
# Load Libraries and Functions #
require(corrplot) ; require(ape)
library(geigen)
source("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/Baypass_executables/baypass_utils.R")
library(tidyverse)
library(GenWin)
library(R.utils)
```
The correlation matrices calculated with the different datasets will be compared to make sure of the consistency in population history reconstruction
```{R}
# Iterate for every combination of 2 numbers between 1 and 50 (pairwise dataset comparison)
omegas <- data.frame()
for (n in 1:length(combn(50,2))){

  i <- combn(50,2)[,n][1]
  k <- combn(50,2)[,n][2]

  assign("omega1",
         as.matrix(read.table(paste0("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/CORE_", i, "_mat_omega.out"))))
  assign("omega2",
         as.matrix(read.table(paste0("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/CORE_", k, "_mat_omega.out"))))

  # Calculate FMD distance value between matrixes
  fmd <- fmd.dist(omega1,omega2)

  # add dataframe entry
  omegas <- rbind(omegas, data.frame(O1 = i, O2 = k, FMD = fmd))
}
summary(omegas)
```
Since the FMD values are all low, all the correlation matrices are very similar between each other. We can take the first as an example and represent it grafically
```{R}
## REPRESENTATION OF CORRELATION MATRICES
# Upload estimate of omega1 and omega2 as trial #
omega1=as.matrix(read.table("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/CORE_1_mat_omega.out"))
omega2=as.matrix(read.table("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/CORE_2_mat_omega.out"))

pop.names=c("CA","KI","LA","MO","TU","UR","VL","YA")
dimnames(omega1)=list(pop.names,pop.names)
dimnames(omega2)=list(pop.names,pop.names)

# Compute and visualize the correlation matrix #
cor.mat1=cov2cor(omega1)
corrplot(cor.mat1,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))
cor.mat2=cov2cor(omega2)
corrplot(cor.mat1,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))

# Visualize the correlation matrix as hierarchical clustering tree #
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat1**2)))
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat2**2)))
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
```
In order to calculate the significativity threshold of the XtX values, I will simulate a dataset of 10k SNPs using one of my Omega correlation matrix. Running BayPass on this neutral simulated dataset will allow me to identify a 99% threshold.
```{R}
# get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution
pi.beta.coef=read.table("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/CORE_1_summary_beta_params.out",h=T)$Mean
# upload the original data to obtain total allele count
allele.data <- geno2YN("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/all.allelecounts")
# Create the POD
simu.all<-simulate.baypass(omega.mat=omega1,nsnp=10000,sample.size=allele.data$NN,
                           beta.pi=pi.beta.coef,pi.maf=0,suffix="PODS")
```
This generates 3 files, one of which (G.PODS) will be used to re-run BayPass and get the XtX values for the simulated dataset:
```{bash}
# Copy the allele file to the EBD server
scp /Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/G.PODS \
ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/BayPass/AlleleCounts/

g_baypass -npop 8 -gfile /home/ebazzicalupo/BayPass/AlleleCounts/G.PODS \
-outprefix /home/ebazzicalupo/BayPass/OutPut/PODS -nthreads 20

scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/BayPass/OutPut/PODS_* \
/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/
```
Now with the re-run analysis we can define a 99% threshold for XtX outliers
```{R}
# get the pod XtX
pod.xtx=read.table("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/PODS_summary_pi_xtx.out",h=T)$M_XtX
# compute the >99.9% threshold
pod.thresh=as.numeric(quantile(pod.xtx,probs=0.999,na.rm=T))
```
The XtX values calculated from all of the datasets can be plotted like this:
```{R}
## OUTLIER SNP DETECTION
# Estimates of the XtX differentiation measures #

# Create a joint table of all the datasets.
# To know which SNP is which I will give each a number based on the actual SNP number of the general dataset:
# dataset1: SNP-1, SNP-51, SNP-101...; dataset2: SNP-2, SNP-52, SNP-102, ...; etc.
# create total snps table
snps.table <- data.frame()

for (n in 1:50){
 # upload the dataset table:
 anacore.snp.res=read.table(paste0("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/CORE_",n,"_summary_pi_xtx.out"),h=T)

 # calculate the database SNP number sequence - 2100553 is the total number of SNPs in all datasets
 anacore.snp.res <- data.frame(anacore.snp.res, SNPnum = seq(n, 2100553, by = 50))

 # add the dataset rows to the total snps table
 snps.table <- rbind(snps.table, anacore.snp.res)
}

# order the table based on the SNP number
snps.table <- snps.table %>% arrange(SNPnum)

# Add SNP ID information
SNPIDs <- read_tsv("/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/BayPass_OutPut/finalset.maf5pc.SNPIDs",
                   col_names = F)[,2-3] %>%
    rename("scaffold" =  X2, "position" = X3)

# Add SNP IDs to the total snps table
snps.table <- data.frame(snps.table,SNPIDs)

# Add Odd/Even column for Manhattan plot
CHR <- data.frame()
for (n in 1:length(unique(snps.table$scaffold))){
   scaffold_lines <- snps.table %>% filter(scaffold==unique(snps.table$scaffold)[n])
   if((n %% 2) == 0) {
    num <- "Even"
   } else {
    num <- "Odd"
   }
   number <- data.frame(colora = rep(num, nrow(scaffold_lines)))
   CHR <- data.frame(rbind(CHR, number))
}
snps.table <- data.frame(cbind(snps.table, CHR))

# Save the table
write.table(x = snps.table,
            file = paste0("BayPass_results/CORE_XtX_results.tsv"),
            quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

# plot XtX values and add the threshold calculated with PODS
manhplot <- ggplot(snps.table, aes(x = SNPnum, y = M_XtX, color = colora)) +
  geom_point(alpha = 0.75, stat = "identity", size=1) +
  scale_color_manual(values= c("Black","darkgray")) +
  geom_hline(yintercept=pod.thresh,linetype="dashed", size=0.5, color="red")
manhplot

# plot XtX values above 15 (pruned) and add the threshold calculated with PODS
pruned.snps.table <- subset(snps.table, M_XtX > 15)
pruned.manhplot <- ggplot(pruned.snps.table, aes(x = SNPnum, y = M_XtX,
                                 color = as.factor(scaffold))) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept=pod.thresh,linetype="dashed", size=0.5)
print(pruned.manhplot)

# try with representing just one scaffold with physical distance between SNPs
sc11arrow.snps.table <- subset(snps.table, scaffold == "scaffold_11_arrow_ctg1")
sc11arrow.manhplot <- ggplot(sc11arrow.snps.table, aes(x = position, y = M_XtX,
                                 color = as.factor(scaffold))) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept=pod.thresh,linetype="dashed", size=0.5)
print(sc11arrow.manhplot)
```
To get the list of SNPs that exceed the threshold value:
```{R}
snps.outliers <- subset(snps.table, M_XtX > pod.thresh)
nrow(snps.outliers)
```

## AUX Model

Now I can run BayPass under the The Auxiliary Variable Covariate model (AUX), in order to check for SNPs associated with a particular environmental co-variable (see Climatic_Variables.md). First I need to divide my table with all the data into tables each containing data for one covariable only:
```{bash}
cd /home/ebazzicalupo/BayPass
varLIST=($(cut -f1 Covariate_Data/fivepcs_populations_data_matrix.tsv | grep -v "variable"))

for var in ${varLIST[@]}
 do
  echo "creating ${var} table"
  grep -w "${var}" Covariate_Data/fivepcs_populations_data_matrix.tsv | cut -f2-  | tr '\t' ' ' \
  > Covariate_Data/${var}_data.txt
done
```
The analysis has to be run for each co-variable separately, also I will divide my SNP dataset into 50 as I did before.

After running the first variables as a trial, I have an alternative version of the script (baypass_aux_v2.sh):
```{bash}
cat baypass_aux_v2.sh
# #!/bin/bash
# 
# WD=/home/ebazzicalupo/BayPass
# 
# for n in {1..50}
#  do
#   echo "dataset ${n}"
#   g_baypass -npop 8 -gfile $WD/AlleleCounts/all.allelecounts.${n} -efile $WD/Covariate_Data/${1}_data.txt \
#   -auxmodel -nthreads 3 -scalecov -omegafile $WD/OutPut/CORE_1_mat_omega.out -outprefix $WD/OutPut/AUX_${1}_${n}
# done
```
which will be launched for each variable (PC), running the analysis on the 50 datasets (long time but you can forget about it and wait it to end).
```{bash}
cd /home/ebazzicalupo/BayPass
varLIST=($(cut -f1 Covariate_Data/fivepcs_populations_data_matrix.tsv | grep -v "variable"))
for var in ${varLIST[@]}
 do
  echo "analyzing association with ${var}"
  screen -dmS aux_${var}  sh -c "/home/ebazzicalupo/BayPass/baypass_aux_v2.sh ${var}; exec /bin/bash"
done
```
Download results to laptop to visualize and analyze in R
```{bash}
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/BayPass/OutPut/AUX_PC\* ~/Documents/Selection_Eurasian_Lynx_v2/3-Identify_Candidate_Loci/tables/BayPass_OutPut_pcs
```
To plot the results in R:
```{R}
library(tidyverse)
variables <- c("PC1", "PC2", "PC3", "PC4", "PC5")

for (i in 1:length(variables)){

 var <- variables[i]

 snps.table <- data.frame()

 for (n in 1:50){
  # upload the dataset table:
  var.snp=read.table(paste0("~/Documents/Selection_Eurasian_Lynx_v2/3-Identify_Candidate_Loci/tables/BayPass_OutPut_pcs/AUX_",
                            var,"_",n,"_summary_betai.out"),h=T)

  # calculate the database SNP number sequence - 2100553 is the total number of SNPs in all datasets
  var.snp <- data.frame(var.snp, SNPnum = seq(n, 2100553, by = 50))

  # add the dataset rows to the total snps table
  snps.table <- rbind(snps.table, var.snp)
 }

 # order the table based on the SNP number
 snps.table <- snps.table %>% arrange(SNPnum)

 # Add SNP ID information
 SNPIDs <- read_tsv("3-Identify_Candidate_Loci/tables/BayPass_OutPut/finalset.maf5pc.SNPIDs", col_names = F)[,2-3] %>%
    rename("scaffold" =  X2, "position" = X3)

 # Add SNP IDs to the total snps table
 snps.table <- data.frame(snps.table,SNPIDs)

 # Add chromosome numbers and Odd/Even for colors
 CHRnumber <- data.frame()
 for (n in 1:length(unique(snps.table$scaffold))){
   scaffold_lines <- snps.table %>% filter(scaffold==unique(snps.table$scaffold)[n])
   if((n %% 2) == 0) {
    num <- "Even"
   } else {
    num <- "Odd"
   }
   number <- data.frame(chromosome = rep(n, nrow(scaffold_lines)), colora = num)
   CHRnumber <- data.frame(rbind(CHRnumber, number))
 }
 snps.table <- data.frame(cbind(snps.table, CHRnumber))

 # Plot BF - ADJUST LEGEND AND X axis (show chromosome/scaffold)
 manhplot <- ggplot(snps.table, aes(x = SNPnum, y = BF.dB., color = colora)) +
   geom_point(alpha = 0.75, stat = "identity") +
   scale_color_manual(values= c("Black","Grey")) +
   theme_light()
    # TITLE
    # LEGEND
    # X AXIS SCAFFOLD NAMES BELOW
   # geom_hline(yintercept=20,linetype="dashed", size=0.5, color="red")

 # save plot
 ggsave(filename = paste0("3-Identify_Candidate_Loci/plots/",var,"_BayPass_manhattanplot.pdf"),
        plot=manhplot, height=8, width = 12, units ="cm", dpi="print")

 # Outliers Table
 # snps.outliers <- subset(snps.table, BF.dB. > 20) %>%
 #   select(scaffold, position, SNPnum, BF.dB.)
 # write.table(x = snps.outliers,file = paste0("BayPass_results/",var,"_outliers_SNPs.tsv"),quote=FALSE,
 #             col.names = T, row.names = FALSE, sep= "\t")
}
```
