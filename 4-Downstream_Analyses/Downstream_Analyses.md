---
title: "Downstream_Analyses"
author: "Enrico"
date: "7/12/2021"
output: html_document
editor_options:
  chunk_output_type: console
---

Downstream analyses of the candidate loci include different approaches to explore mainly 3 different aspects:
  
  1. Spatial gradients of environmental effects on genetics
  2. Effects of selection on the genetic structure of populations
  3. Functional annotation of the candidate loci

## 1. Spatial gradients of environmental effects on genetics

Generalized Dissimilarity Modeling (GDM) will be run to analyze how environmental pressure gradients are affecting genetics differentially in space.

These analyses are based on a mix of these two working examples:
https://github.com/pgugger/LandscapeGenomics/blob/master/2017/Exercise4.md#generalized-dissimilarity-modeling-gdm
https://www.mountainmanmaier.com/software/pop_genom/#generalized-dissimilarity-modeling

To get the SNPs in a format readable in R (not plink's RAW):
```{bash}
# on genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/
# get VCF of topsnps only
bedtools intersect -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
 -b Intersect/total_intersect_candidate_windows_topsnps.bed \
 > VCF/ll_allsamples_total_intersect_candidate_windows_topsnps.vcf

cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF

# get matrix with rows=samples, cols=SNPs, values=number of alternative alleles in sample(0,1,2)
vcftools --vcf ll_allsamples_total_intersect_candidate_windows_topsnps.vcf --012 \
 --out topsnps_rformat

# edit outputs to create a table
cut -f2- topsnps_rformat.012 | sed 's/-1/NA/g' >topsnps_rformat.temp
tr -d '\t' <topsnps_rformat.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
paste <(echo "ID" | cat - topsnps_rformat.012.indv) <(echo "" | cat header - topsnps_rformat.temp) > topsnps_rformat.forR
rm header topsnps_rformat.temp
```
Download to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/topsnps_rformat.forR Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
Prepare R
```{R}
# load libraries
library(tidyverse)
library(gdm)
```
Import data
```{R}
# Samples to be excluded
samples_to_remove <- c("c_ll_ba_0216","c_ll_ba_0233","c_ll_cr_0211","h_ll_ba_0214","h_ll_ba_0215")

# SNP data
snp <- read.table("4-Downstream_Analyses/tables/topsnps_rformat.forR", header = T, row.names = 1)
snp <- snp[!(row.names(snp) %in% samples_to_remove),]

# env.data
clim.data <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")

clim.data <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")

# coordinates
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_wholeset.csv",
                          col_names = T, delim = ';')  %>% column_to_rownames(., var="id")
coord_table <- coord_table[!(row.names(coord_table) %in% samples_to_remove),]

# dataframe with data + coordinates
clim.points <- data.frame(clim.data, x=as.numeric(coord_table$longitude), y=as.numeric(coord_table$latitude)) %>% 
  rownames_to_column(., var="ID")
```

```{R}
snp.dist <- dist(snp, diag = T, upper = T)   # generate Euclidean distance matrix
snp.dist.1 <- snp.dist/(max(snp.dist))  #rescale by dividing by max value
snp.dist.1 <- as.matrix(snp.dist.1)
snp.dist.1 <- cbind(ID = rownames(snp.dist.1), snp.dist.1)  # make the row names an actual column named "ID"
rownames(snp.dist.1) <- NULL  #remove prior row names
```

```{R}
gdm.input <- gdm::formatsitepair(bioData=snp.dist.1, bioFormat=3, 
                            predData=clim.points, siteColumn="ID", 
                            XColumn="x", YColumn="y")
gdm <- gdm(gdm.input, geo = T, splines = NULL, knots = NULL)

summary(gdm)
gdm$explained

gdm.importance <- gdm.varImp(gdm.input, geo=T, nPerm=100, parallel=TRUE, cores=32)

pdf(file = paste0("4-Downstream_Analyses/plots/GDM_importance.pdf"),
    width = 8,
    height = 8)
barplot(sort(gdm.importance[[2]][,1], decreasing=T))
dev.off()

#               fullModel  fullModel-1  fullModel-2 fullModel-3 fullModel-4 fullModel-5 fullModel-6
# Geographic 1.887268e+01 1.947259e+01 1.947259e+01 23.76652482   23.796720   36.262058    60.89688
# bio2       6.985516e+00 6.985516e+00 6.985516e+00  6.98463942    6.926600    6.635340          NA
# bio5       8.864149e-02 8.864149e-02 8.864149e-02  0.09357642          NA          NA          NA
# bio6       1.200511e+00 1.200511e+00 1.200511e+00  1.19958041    1.500752          NA          NA
# bio8       0.000000e+00 0.000000e+00           NA          NA          NA          NA          NA
# bio13      9.419809e-04 9.419809e-04 9.419809e-04          NA          NA          NA          NA
# jan_depth  0.000000e+00           NA           NA          NA          NA          NA          NA
# snow_days  5.106221e+00 6.140592e+00 6.140592e+00  6.13970826    6.080373    8.598162    11.44324


plot(gdm, plot.layout = c(2, 5))

pdf(file = paste0("4-Downstream_Analyses/plots/GDM_effect_uncertainty.pdf"),
    width = 8,
    height = 8)
plotUncertainty(gdm.input, sampleSites=0.70, bsIters=100, geo=T, plot.layout=c(2,4))
dev.off()
```
Transform climate data layers based on modeled results:
```{R}
library(raster)
library(rgdal)
library(grDevices)
library(sf)

# bio variables
worldclim <- getData("worldclim", var = "bio", res = 10)
DEM2 <- worldclim[[2]]
DEM5 <- worldclim[[5]]
DEM6 <- worldclim[[6]]
DEM8 <- worldclim[[8]]
DEM13 <- worldclim[[13]]

# january mean depth
jan_mean_depth_stack <- raster()
for (y in 1999:2018){
  print(y)
  data <- stack(paste0("2-Prepare_Environmental_Data/tables/Snow_data_daily/cmc_sdepth_dly_", y, "_v01.2.tif"))
  data_jan <- data[[c(1:31)]]
  mean_jan <- calc(data_jan, fun = mean)
  jan_mean_depth_stack <- stack(jan_mean_depth_stack, mean_jan)
}
jan_mean_depth <- calc(jan_mean_depth_stack, fun = mean)
plot(jan_mean_depth)
writeRaster(jan_mean_depth, filename = "2-Prepare_Environmental_Data/tables/jan_mean_depth.tif")

jan_mean_depth <- raster("2-Prepare_Environmental_Data/tables/jan_mean_depth.tif")

jan_mean_depth_new <- projectRaster(jan_mean_depth, worldclim)

# mean number of yearly snow days 
snow_days_stack <- raster()
for (y in 1999:2018){
  print(y)
  data <- stack(paste0("2-Prepare_Environmental_Data/tables/Snow_data_daily/cmc_sdepth_dly_", y, "_v01.2.tif"))
  binary_data <- data > 0
  sum_snow <- calc(data_jan, fun = sum)
  snow_days_stack <- stack(snow_days_stack, sum_snow)
}
mean_snow_days <- calc(snow_days_stack, fun = mean)
plot(mean_snow_days)
writeRaster(mean_snow_days, filename = "2-Prepare_Environmental_Data/tables/mean_snow_days.tif")

mean_snow_days <- raster("2-Prepare_Environmental_Data/tables/mean_snow_days.tif")

mean_snow_new <- projectRaster(mean_snow_days, worldclim)

# combine layers
clim.layer <- stack(DEM2,DEM5,DEM6,DEM8,DEM13,jan_mean_depth_new,mean_snow_new)

extent <- c(-10, 180, 20, 90) 
clim.layer.crop <- crop(clim.layer, extent)

clim.trans <- gdm.transform(gdm, clim.layer.crop)
clim.rast <- na.omit(getValues(clim.trans))

pdf("4-Downstream_Analyses/plots/GDM_vars_Maps.pdf")
plot(clim.trans)
dev.off()

plot(clim.trans[[7]])
points(clim.points$x,clim.points$y, pch=4, cex=0.2)
plot(clim.trans[[3]])
points(clim.points$x,clim.points$y, pch=4, cex=0.2)

pca <- prcomp(clim.rast)
pca.rast <- predict(clim.trans, pca, index=1:3)
pca.rast[[1]] <- (pca.rast[[1]]-pca.rast[[1]]@data@min) / (pca.rast[[1]]@data@max-pca.rast[[1]]@data@min)*255
pca.rast[[2]] <- (pca.rast[[2]]-pca.rast[[2]]@data@min) / (pca.rast[[2]]@data@max-pca.rast[[2]]@data@min)*255
pca.rast[[3]] <- (pca.rast[[3]]-pca.rast[[3]]@data@min) / (pca.rast[[3]]@data@max-pca.rast[[3]]@data@min)*255

plotRGB(pca.rast, r=1, g=2, b=3, bgalpha=0)

plot(pca.rast[[3]])

balkans <- c(10, 35, 30, 50) 
pca.rast.balkans <- crop(pca.rast, balkans)

plotRGB(pca.rast.balkans, r=1, g=2, b=3, bgalpha=0)


points(clim.points$x,clim.points$y, pch=4, cex=0.2)

```

## 2. Effects of selection on the genetic structure of populations

To explore this aspect we will run RDA again, but this time we will compare structure recovered from selected vs unselected loci. We will be investigating only the variables that were modeled to have a strong effect on genetic data by GDM (Bio2 and Snow Days).

To get the TOPSNP of each candidate window of the strong effect variables:
```{bash}
# on genomics-a:
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/

# get topsnps from all uncorrelated predictors
rm Intersect/topvars_topsnps.range
for var in bio2 snow_days
 do
  echo "${var}"
  cat GenWin/${var}_topsnps.range >> Intersect/topvars_topsnps.range
done
```
Problem = more than one topsnp for each candidate window (window is same topsnp is different because from different variable or originally window was split in 2)
Solution = get topsnps from each candidate window with only one topnsp and then get one random topsnp from each window with more than one
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/Intersect

# get list of candidate windows with ONE topsnp
bedtools intersect -wo -a total_intersect_candidate_windows.bed \
 -b <(cat topvars_topsnps.range | sort -k 1,1 -k2,2n) |
 cut -f1,2,3 | uniq -u > topvars_topsnps_nodupwindows.bed
 
# get topsnp from each window with ONE topsnp
bedtools intersect -a <(cat topvars_topsnps.range | sort -k 1,1 -k2,2n) \
 -b topsnps_nodupwindows.bed > topvars_topsnps_nodupwindows_onesnp.bed

# get list of candidate windows with MORE than one topsnp
bedtools intersect -wo -a total_intersect_candidate_windows.bed \
 -b <(cat topvars_topsnps.range | sort -k 1,1 -k2,2n) |
 cut -f1,2,3 | uniq -d > topvars_topsnps_dupwindows.bed

# get only one topsnp for each of the windows with duplicates
while read p; do
  bedtools intersect -a <(cat topvars_topsnps.range | sort -k 1,1 -k2,2n) -b <(echo "$p") | shuf -n 1
done < topvars_topsnps_dupwindows.bed > topvars_topsnps_dupwindows_onerand.bed

# get full list of topsnps (only one per window)
cat topvars_topsnps_nodupwindows_onesnp.bed topvars_topsnps_dupwindows_onerand.bed | 
 sort -k 1,1 -k2,2n \
 > topvars_intersect_candidate_windows_topsnps.bed
```

Extract VCF and PLINK of the topsnps
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/

bedtools intersect -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
 -b Intersect/topvars_intersect_candidate_windows_topsnps.bed \
 > VCF/ll_allsamples_topvars_intersect_candidate_windows_topsnps.vcf

cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/
# I manually created a file listing samples to remove for plink with following format:
# c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

plink_1.9 --vcf ll_allsamples_topvars_intersect_candidate_windows_topsnps.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove samplestoremove.txt --geno 0 \
 --recode A --out topvars_intersect_candidate_windows_topsnps
```
Copy to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/topvars_intersect_candidate_windows_topsnps.raw Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
Prepare R for RDA on laptop
```{R}
# load libraries
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(adegenet)

# Load environmental data
env.predictors <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")

env.predictors <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")
env.predictors <- env.predictors[, c('bio2', 'snow_days')]

# Load GenoType Data - choose set of vars file
gt_data <- read.PLINK("4-Downstream_Analyses/tables/topvars_intersect_candidate_windows_topsnps.raw")
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
pdf(file = paste0("4-Downstream_Analyses/plots/candidate_rda1_rda2.pdf"),
    width = 8,
    height = 8)

plot(rda, type="n", scaling=3)
#points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num]) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1)     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

dev.off()

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

## RDA on neutral regions only

Filter the VCF to remove candidate regions and all genes.

Intergenic regions were identified by Dani as follows:
```{bash}
cd /GRUPOS/grupolince/reference_genomes/lynx_canadensis

awk -F"\t" '$3 == "gene" {printf ("%s\t%s\t%s\n", $1, $4-5001, $5+5000)}' lc4.NCBI.nr_main.gff3 \
 > lc4.NCBI.nr_main.genes.plus5000.temp_bed

join -1 1 -2 1 <(LANG=en_EN sort -k1,1 -k2,2n -k3,3n lc_ref_all_the_genome.bed) \
 <(LANG=en_EN sort -k1,1 -k2,2n -k3,3n lc4.NCBI.nr_main.genes.plus5000.temp_bed) | 
 awk -v OFS='\t' '{if ($4<0) {$4="1"} else {$4=$4}; print}' | 
 awk -v OFS='\t' '{if ($5>$3) {$5=$3} else {$5=$5}; print $1,$4,$5}' | 
 bedtools merge -i stdin -d 1 > lc4.NCBI.nr_main.genes.plus5000.bed

rm lc4.NCBI.nr_main.genes.plus5000.temp_bed

bedtools subtract -a lc_ref_all_the_genome.bed -b lc4.NCBI.nr_main.genes.plus5000.bed \
 > lc4.NCBI.nr_main.intergenic.buffer5000.bed
```
Remove those regions from the VCF
```{bash}
# on genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF

bedtools subtract -header -a ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
 -b /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.intergenic.buffer5000.bed \
 > ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.vcf
``` 
From 4983054 whole genome SNPs to 2507494 intergenic SNPs.

Now remove also any candidate of selection - the SNPs from the RDA and the GenWin windows from BayPass
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/
# create a new VCF to substract from sequentially
cp VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.vcf \
 VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf
 
# Subtract GenWin windows sequentially from new VCF
for var in bio2 bio5 bio6 bio8 bio13 jan_depth snow_days
 do
  echo "filtering ${var} regions"
  bedtools subtract -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
   -b GenWin/${var}_GenWin_windows_outliers.bed \
   > tmp && mv tmp VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf
done

# Subtract RDA SNPs from new VCF
bedtools subtract -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
 -b RDA/rda_candidate_snps.bed \
 > tmp && mv tmp VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf
```
We have a total of 2284573 non-selected intergenic SNPs

Now to create a RAW file with plink with only 5k uncorrelated intergenic non-selected SNPs
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/
# I manually created a file listing samples to remove for plink with following format:
# c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

# prune snps based on ld (VIF < 2)
plink_1.9 --vcf ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# --geno 0 \
--remove samplestoremove.txt --indep 100 10 2

# take only a subsample of 10k of the pruned SNPs
plink_1.9 --vcf ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--remove samplestoremove.txt --geno 0 --extract <(shuf -n 10000 plink.prune.in) \
--recode A --out intergenic_neutral_snps
```
Copy to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/intergenic_neutral_snps.raw Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
Prepare R for RDA on laptop
```{R}
# load libraries
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(adegenet)

# Load environmental data
env.predictors <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")
env.predictors <- env.predictors[, c('bio2', 'snow_days')]

# Load GenoType Data - choose set of vars file
gt_data <- read.PLINK("4-Downstream_Analyses/tables/intergenic_neutral_snps.raw")
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
pdf(file = paste0("4-Downstream_Analyses/plots/neutral_rda1_rda2.pdf"),
    width = 8,
    height = 8)

plot(rda, type="n", scaling=3)
points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num]) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1)     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

dev.off()

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

## Discounting geographic distance from RDA

```{R}
library(geosphere)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(adegenet)

### build geographic distance matrix

samples_to_remove <- c("c_ll_ba_0216","c_ll_ba_0233","c_ll_cr_0211","h_ll_ba_0214","h_ll_ba_0215")

coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_wholeset.csv",
                          col_names = T, delim = ';')  %>% column_to_rownames(., var="id")
coord_table <- coord_table[!(row.names(coord_table) %in% samples_to_remove),]

coords <- data.frame(x=as.numeric(coord_table$longitude), y=as.numeric(coord_table$latitude),
                     row.names = rownames(coord_table))

distm <- distm(as.matrix(coords), fun = distVincentyEllipsoid)
d <- vegdist(distm)
### RDA discounting variance explained by geo distance matrix

# Load environmental data
env.predictors <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")

env.predictors <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")
env.predictors <- env.predictors[, c('bio2', 'snow_days')]

env.predictors <- data.frame(bio2=env.predictors$bio2)

env.dist <- vegdist(env.predictors, upper=T, diag=T)

# Load GenoType Data - choose set of vars file
gt_data <- read.PLINK("4-Downstream_Analyses/tables/topvars_intersect_candidate_windows_topsnps.raw")
gt_data_tsv <- data.frame(as.matrix(gt_data))

gen.dist <- vegdist(gt_data_tsv)

# Run RDA
rda <- rda(gt_data ~ . + Condition(coords$x) + Condition(coords$y), data=env.predictors, scale=T)
signif.axis <- anova.cca(rda, by="axis", parallel=getOption("mc.cores"))


rda <- rda(gt_data, env.predictors, coords)
rda <- rda(gt_data, coords)
rda <- rda(gt_data ~ ., data=coords)

biplot(rda, scaling = 0)

# DB RDA

dbrda <- dbrda(as.matrix(gen.dist) ~ as.matrix(env.dist) + Condition(as.matrix(d)), scale=T)
screeplot(dbrda)
summary(eigenvals(dbrda, model = "constrained"))
```

## PCA of candidates

Run PCA with plink on genomics-a
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/
# I manually created a file listing samples to remove for plink with following format:
# c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

plink_1.9 --vcf ll_allsamples_topvars_intersect_candidate_windows_topsnps.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove samplestoremove.txt \
 --pca --out topvars_intersect_candidate_windows_topsnps_pca
```
Copy eigen files to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/topvars_intersect_candidate_windows_topsnps_pca.eigen\* Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
prepare R
```{R}
library(tidyverse)
library(RColorBrewer)
library (viridis)
```
plot PCA
```{R}
# import eigen vec and val
pca <- read_table2("4-Downstream_Analyses/tables/topvars_intersect_candidate_windows_topsnps_pca.eigenvec",
                   col_names = FALSE)
eigenval <- scan("4-Downstream_Analyses/tables/topvars_intersect_candidate_windows_topsnps_pca.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# add population
loc <- rep(NA, length(pca$ind))
loc[grep("ca", pca$ind)] <- "#B8860b"
loc[grep("po", pca$ind)] <- viridis_pal()(5)[3]
loc[grep("ba", pca$ind)] <- "#A035AF"
loc[grep("cr", pca$ind)] <- brewer.pal(12,"Paired")[9]
loc[grep("no", pca$ind)] <- viridis_pal()(5)[2]
loc[grep("ka", pca$ind)] <- brewer.pal(12,"Paired")[7]
loc[grep("ki", pca$ind)] <- viridis_pal()(5)[1]
loc[grep("la", pca$ind)] <- brewer.pal(12,"Paired")[3]
loc[grep("to", pca$ind)] <- brewer.pal(12,"Paired")[7]
loc[grep("og", pca$ind)] <- brewer.pal(12,"Paired")[7]
loc[grep("tu", pca$ind)] <- brewer.pal(12,"Paired")[8]
loc[grep("ur", pca$ind)] <- "#0F4909"
loc[grep("vl", pca$ind)] <- brewer.pal(12,"Paired")[5]
loc[grep("ya", pca$ind)] <- brewer.pal(12,"Paired")[6]

# remake data.frame
pca <- as.tibble(data.frame(pca, loc))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# then make a plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# plot PCA without labels
pdf(file = paste0("4-Downstream_Analyses/plots/candidate_pca1_pca2.pdf"),
    width = 8,
    height = 8)

ggplot() + 
  geom_point(data=pca, aes(PC1, PC2), size = 4, pch=21, bg=loc, col="gray32") +
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  scale_color_manual(values=cols) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

dev.off()
```

## PCA of neutral regions

Run PCA with plink on genomics-a
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/
# I manually created a file listing samples to remove for plink with following format:
# c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9
# c_ll_ca_0245 c_ll_ca_0245 0 0 0 -9
# c_ll_ca_0248 c_ll_ca_0248 0 0 0 -9
# c_ll_ca_0254 c_ll_ca_0254 0 0 0 -9
# c_ll_ca_0247 c_ll_ca_0247 0 0 0 -9
# c_ll_ca_0259 c_ll_ca_0259 0 0 0 -9

plink_1.9 --vcf ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove samplestoremove_pca.txt --extract <(shuf -n 5000 plink.prune.in) \
 --pca --out intergenic_neutral_snps_pca
```
Copy eigen files to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/intergenic_neutral_snps_pca.eigen\* Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
prepare R
```{R}
library(tidyverse)
library(RColorBrewer)
library (viridis)
```
plot PCA
```{R}
# import eigen vec and val
pca <- read_table("4-Downstream_Analyses/tables/intergenic_neutral_snps_pca.eigenvec",
                   col_names = FALSE)
eigenval <- scan("4-Downstream_Analyses/tables/intergenic_neutral_snps_pca.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# add population
loc <- rep(NA, length(pca$ind))
loc[grep("ca", pca$ind)] <- "#B8860b"
loc[grep("po", pca$ind)] <- viridis_pal()(5)[3]
loc[grep("ba", pca$ind)] <- "#A035AF"
loc[grep("cr", pca$ind)] <- brewer.pal(12,"Paired")[9]
loc[grep("no", pca$ind)] <- viridis_pal()(5)[2]
loc[grep("ka", pca$ind)] <- brewer.pal(12,"Paired")[7]
loc[grep("ki", pca$ind)] <- viridis_pal()(5)[1]
loc[grep("la", pca$ind)] <- brewer.pal(12,"Paired")[3]
loc[grep("to", pca$ind)] <- brewer.pal(12,"Paired")[7]
loc[grep("og", pca$ind)] <- brewer.pal(12,"Paired")[7]
loc[grep("tu", pca$ind)] <- brewer.pal(12,"Paired")[8]
loc[grep("ur", pca$ind)] <- "#0F4909"
loc[grep("vl", pca$ind)] <- brewer.pal(12,"Paired")[5]
loc[grep("ya", pca$ind)] <- brewer.pal(12,"Paired")[6]

# remake data.frame
pca <- as.tibble(data.frame(pca, loc))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# then make a plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# plot PCA without labels
pdf(file = paste0("4-Downstream_Analyses/plots/neutral_pca1_pca2.pdf"),
    width = 8,
    height = 8)

ggplot() + 
  geom_point(data=pca, aes(PC1, PC2), size = 4, pch=21, bg=loc, col="gray32") +
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  scale_color_manual(values=cols) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

dev.off()
```
