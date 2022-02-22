---
title: "Downstream_Analyses_pcs"
author: "Enrico"
date: "2/21/2022"
output: html_document
editor_options:
  chunk_output_type: console
---
Downstream analyses of the candidate loci include different approaches to explore mainly 3 different aspects:

1.  Spatial gradients of environmental effects on genetics
2.  Effects of selection on the genetic structure of populations
3.  Functional annotation of the candidate loci

## 1. Spatial gradients of environmental effects on genetics

Generalized Dissimilarity Modeling (GDM) will be run to analyze how environmental pressure gradients are affecting genetics differentially in space.

These analyses are based on a mix of these two working examples: <https://github.com/pgugger/LandscapeGenomics/blob/master/2017/Exercise4.md#generalized-dissimilarity-modeling-gdm> <https://www.mountainmanmaier.com/software/pop_genom/#generalized-dissimilarity-modeling>

To get the SNPs in a format readable in R (not plink's RAW):

```{bash}
# on genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/
# get VCF of topsnps only
bedtools intersect -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
 -b Intersect/total_intersect_candidate_fivepcs_windows_topsnps.bed \
 > VCF/ll_allsamples_total_intersect_candidate_fivepcs_windows_topsnps.vcf

cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF

# get matrix with rows=samples, cols=SNPs, values=number of alternative alleles in sample(0,1,2)
vcftools --vcf ll_allsamples_total_intersect_candidate_fivepcs_windows_topsnps.vcf --012 \
 --out fivepcs_topsnps_rformat

# edit outputs to create a table
cut -f2- fivepcs_topsnps_rformat.012 | sed 's/-1/NA/g' >fivepcs_topsnps_rformat.temp
tr -d '\t' <fivepcs_topsnps_rformat.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
paste <(echo "ID" | cat - fivepcs_topsnps_rformat.012.indv) <(echo "" | cat header - fivepcs_topsnps_rformat.temp) \
 > fivepcs_topsnps_rformat.forR
rm header fivepcs_topsnps_rformat.temp
```
Download to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/fivepcs_topsnps_rformat.forR Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
Prepare R
```{R}
# load libraries
library(tidyverse)
library(gdm)
library(raster)
library(rgdal)
library(grDevices)
library(sf)
library(viridis)
library(RColorBrewer)
```
Import data
```{R}
# Samples to be excluded
samples_to_remove <- c("c_ll_ba_0216","c_ll_ba_0233","c_ll_cr_0211","h_ll_ba_0214","h_ll_ba_0215")

# SNP data
snp <- read.table("4-Downstream_Analyses/tables/fivepcs_topsnps_rformat.forR", header = T, row.names = 1)
snp <- snp[!(row.names(snp) %in% samples_to_remove),]

# env.data all variables
varmatrix <- read_tsv("2-Prepare_Environmental_Data/WorldClim_table_persample.tsv") %>% column_to_rownames(., var="sample")
snowdata <- read_tsv("2-Prepare_Environmental_Data/Snow_table_persample.tsv") %>% column_to_rownames(., var="sample")

clim.data <- cbind(varmatrix, snowdata)

# coordinates
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_wholeset.csv",
                          col_names = T, delim = ';')  %>% column_to_rownames(., var="id")
coord_table <- coord_table[!(row.names(coord_table) %in% samples_to_remove),]

# dataframe with data + coordinates
clim.points <- data.frame(clim.data, x=as.numeric(coord_table$longitude), y=as.numeric(coord_table$latitude)) %>% 
  rownames_to_column(., var="ID")
```
Calculate Euclidian distance between samples based on SNP data (dissimilarity)
```{R}
snp.dist <- dist(snp, diag = T, upper = T)   # generate Euclidean distance matrix
snp.dist.1 <- snp.dist/(max(snp.dist))  #rescale by dividing by max value
snp.dist.1 <- as.matrix(snp.dist.1)
snp.dist.1 <- cbind(ID = rownames(snp.dist.1), snp.dist.1)  # make the row names an actual column named "ID"
rownames(snp.dist.1) <- NULL  #remove prior row names
```
Run GDM. Check relative importance of different variables on shaping genetic dissimilarity between samples. 
```{R}
# Prepare input
gdm.input <- gdm::formatsitepair(bioData=snp.dist.1, bioFormat=3, 
                            predData=clim.points, siteColumn="ID", 
                            XColumn="x", YColumn="y")
# Run GDM
gdm <- gdm(gdm.input, geo = T, splines = NULL, knots = NULL)

# Print summary
summary(gdm)
# Overall power of GDM analysis
gdm$explained

# Calculate each variable's importance
gdm.importance <- gdm.varImp(gdm.input, geo=T, nPerm=50, parallel=TRUE, cores=8)

# Create data frame for plot - full model importance values
gdm.imp.df <- data.frame(importance=sort(gdm.importance[[2]][,1], decreasing=T)) %>% rownames_to_column(., var="variable")

gdm.imp.df
# variable importance
#       bio2 9.80252092
#  snow_days 9.54946427
# Geographic 8.09948414
#       bio3 1.42863056
#      bio18 0.63829982
#      bio15 0.30638537
#       bio1 0.18318020
#      bio19 0.01290429
#       bio4 0.00000000
#       bio5 0.00000000
#       bio6 0.00000000
#       bio7 0.00000000
#       bio8 0.00000000
#       bio9 0.00000000
#      bio10 0.00000000
#      bio11 0.00000000
#      bio12 0.00000000
#      bio13 0.00000000
#      bio14 0.00000000
#      bio16 0.00000000
#      bio17 0.00000000
#  jan_depth 0.00000000


# Plot importance in a barplot ordered from most to least important
p <- ggplot(data=gdm.imp.df, aes(x=variable, y=importance)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=c(rev(gdm.imp.df$variable))) +
  theme_minimal() +
  coord_flip()

pdf(file = paste0("4-Downstream_Analyses/plots/GDM_fivepcs_importance.pdf"),
    width = 8,
    height = 8)
p
dev.off()

# Plot the effect function of each variable, with uncertainty in estimation (confidence interval?)
pdf(file = paste0("4-Downstream_Analyses/plots/GDM_fivepcs_effect_uncertainty.pdf"),
    width = 8,
    height = 8)
plotUncertainty(gdm.input, sampleSites=0.70, bsIters=100, geo=T, plot.layout=c(2,4))
dev.off()
```
Transform climate data layers based on modeled results:
```{R}
# load worldclim bio variables
worldclim <- getData("worldclim", var = "bio", res = 10)

# january mean depth raster was created like this:
# jan_mean_depth_stack <- raster()
# for (y in 1999:2018){
#   print(y)
#   data <- stack(paste0("2-Prepare_Environmental_Data/tables/Snow_data_daily/cmc_sdepth_dly_", y, "_v01.2.tif"))
#   data_jan <- data[[c(1:31)]]
#   mean_jan <- calc(data_jan, fun = mean)
#   jan_mean_depth_stack <- stack(jan_mean_depth_stack, mean_jan)
# }
# jan_mean_depth <- calc(jan_mean_depth_stack, fun = mean)
# plot(jan_mean_depth)
# writeRaster(jan_mean_depth, filename = "2-Prepare_Environmental_Data/tables/jan_mean_depth.tif")

# load january depth
jan_mean_depth <- raster("2-Prepare_Environmental_Data/tables/jan_mean_depth.tif")
jan_mean_depth_new <- projectRaster(jan_mean_depth, worldclim)

# mean number of yearly snow days raster was created like this:
# snow_days_stack <- raster()
# for (y in 1999:2018){
#   print(y)
#   data <- stack(paste0("2-Prepare_Environmental_Data/tables/Snow_data_daily/cmc_sdepth_dly_", y, "_v01.2.tif"))
#   binary_data <- data > 0
#   sum_snow <- calc(data_jan, fun = sum)
#   snow_days_stack <- stack(snow_days_stack, sum_snow)
# }
# mean_snow_days <- calc(snow_days_stack, fun = mean)
# plot(mean_snow_days)
# writeRaster(mean_snow_days, filename = "2-Prepare_Environmental_Data/tables/mean_snow_days.tif")

# load mean snow days
mean_snow_days <- raster("2-Prepare_Environmental_Data/tables/mean_snow_days.tif")
mean_snow_new <- projectRaster(mean_snow_days, worldclim)

# combine climate layers
clim.layer <- stack(worldclim,jan_mean_depth_new,mean_snow_new)

# Transform the climate layers based on predicted effects on genetic composition (GDM results)
clim.trans <- gdm.transform(gdm, clim.layer)
```
Load distribution map
```{R}
# distribution borders
distr.map <- readOGR("2-Prepare_Environmental_Data/redlist_species_data/data_0.shp")
```
Add population colors to samples for plotting
```{R}
# Add populations and colors
loc <- rep(NA, NROW(clim.data))
loc[grep("ba", rownames(clim.data))] <- "Balkans"
loc[grep("ca", rownames(clim.data))] <- "Caucasus"
loc[grep("cr", rownames(clim.data))] <- "Carpathians"
loc[grep("ka", rownames(clim.data))] <- "Mongolia"
loc[grep("ki", rownames(clim.data))] <- "Kirov"
loc[grep("la", rownames(clim.data))] <- "Latvia"
loc[grep("no", rownames(clim.data))] <- "Norway"
loc[grep("po", rownames(clim.data))] <- "NE-Poland"
loc[grep("og", rownames(clim.data))] <- "Mongolia"
loc[grep("to", rownames(clim.data))] <- "Mongolia"
loc[grep("tu", rownames(clim.data))] <- "Tuva"
loc[grep("ur", rownames(clim.data))] <- "Urals"
loc[grep("vl", rownames(clim.data))] <- "Vladivostok"
loc[grep("ya", rownames(clim.data))] <- "Yakutia"

num <- rep(NA, NROW(clim.data))
num[grep("ba", rownames(clim.data))] <- 1
num[grep("ca", rownames(clim.data))] <- 3
num[grep("cr", rownames(clim.data))] <- 2
num[grep("ka", rownames(clim.data))] <- 6
num[grep("ki", rownames(clim.data))] <- 4
num[grep("la", rownames(clim.data))] <- 5
num[grep("no", rownames(clim.data))] <- 8
num[grep("po", rownames(clim.data))] <- 7
num[grep("og", rownames(clim.data))] <- 6
num[grep("to", rownames(clim.data))] <- 6
num[grep("tu", rownames(clim.data))] <- 9
num[grep("ur", rownames(clim.data))] <- 10
num[grep("vl", rownames(clim.data))] <- 11
num[grep("ya", rownames(clim.data))] <- 12

ola <- data.frame(sample = rownames(clim.data), pop = loc, n = num)
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
Plot RGB map of 3 most affecting variables (bio2, bio3 and snow_days)
```{R}
# check which layers are which variable:

clim.trans[[4]]@data@names
# bio2

clim.trans[[5]]@data@names
# bio3

clim.trans[[9]]@data@names
# mean_snow_days

# Transform the values in the wanted layers to a scale of 255 for RGB colors
clim.trans[[4]] <- (clim.trans[[4]]-clim.trans[[4]]@data@min) / (clim.trans[[4]]@data@max-clim.trans[[4]]@data@min)*255
clim.trans[[5]] <- (clim.trans[[5]]-clim.trans[[5]]@data@min) / (clim.trans[[5]]@data@max-clim.trans[[5]]@data@min)*255
clim.trans[[9]] <- (clim.trans[[9]]-clim.trans[[9]]@data@min) / (clim.trans[[9]]@data@max-clim.trans[[9]]@data@min)*255

# set the extent of map to plot only interesting part
extent <- c(-10, 180, 20, 90) 
clim.trans.crop <- crop(clim.trans, extent)

# Subset the raster only to include the distribution area (to highlight it)
r2 <- crop(clim.trans, extent(distr.map))
r3 <- mask(r2, distr.map)

# Empty plot to add the graticule of coordinates
empty <- raster(crs=clim.trans)
values(empty) <- 0
empty.crop <- crop(empty, extent)

# Plot bio2 as red, bio3 as green, snow days as blue
pdf("4-Downstream_Analyses/plots/GDM_fivepcs_Rbio2_Gbio3_BsnowD_cool_RGBmap.pdf", width = 14)
plot(empty.crop, col=c('white'), legend=F)
plotRGB(clim.trans.crop, r=4, g=5, b=9, bgalpha=0, alpha=180, add=T)
plotRGB(r3, r=4, g=5, b=9, bgalpha=0, add=T)
plot(distr.map, lwd=1, add=T)
#points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()
```
Plot RGB map of first 3 principal components of a PCA of the raster values
```{R}
# Transform the climate layers based on predicted effects on genetic composition (GDM results)
# RUN AGAIN IF RGB PLOT WAS ALREADY RUN
clim.trans <- gdm.transform(gdm, clim.layer)

# Extract values from the transformed raster
clim.rast <- na.omit(getValues(clim.trans))

# run PCA on raster values
pca <- prcomp(clim.rast)
#                         PC1         PC2           PC3          PC4         PC5
# xCoord          0.993986762  0.05267692  0.0871639792 -0.004127934 -0.01660612
# yCoord         -0.028294824  0.31894731  0.1038583106  0.472037047  0.28856787
# bio1            0.003235923 -0.13516143  0.0329929831 -0.027823340 -0.11920388
# bio2            0.097052778 -0.25383593 -0.9319188118  0.216372070  0.04213606
# bio3           -0.006917332 -0.32360445 -0.0222167738 -0.603895843 -0.15125279
# bio15           0.032905369 -0.09748500 -0.0157423754 -0.425221783  0.66446404
# bio18           0.012305472 -0.02135883  0.0225962729 -0.071441207 -0.55752345
# bio19          -0.008340730  0.01789273  0.0003982097  0.116733104 -0.33837425
# mean_snow_days -0.020208931  0.83534242 -0.3328701335 -0.406618689 -0.10367451

# Create a raster with PC values
pca.rast <- predict(clim.trans, pca, index=1:5)
# Transform the values in the wanted layers to a scale of 255 for RGB colors
pca.rast[[1]] <- (pca.rast[[1]]-pca.rast[[1]]@data@min) / (pca.rast[[1]]@data@max-pca.rast[[1]]@data@min)*255
pca.rast[[2]] <- (pca.rast[[2]]-pca.rast[[2]]@data@min) / (pca.rast[[2]]@data@max-pca.rast[[2]]@data@min)*255
pca.rast[[3]] <- (pca.rast[[3]]-pca.rast[[3]]@data@min) / (pca.rast[[3]]@data@max-pca.rast[[3]]@data@min)*255
pca.rast[[3]] <- (pca.rast[[3]]-pca.rast[[3]]@data@min) / (pca.rast[[3]]@data@max-pca.rast[[3]]@data@min)*255
pca.rast[[5]] <- (pca.rast[[5]]-pca.rast[[5]]@data@min) / (pca.rast[[5]]@data@max-pca.rast[[5]]@data@min)*255

# set the extent of map to plot only interesting part
extent <- c(-10, 180, 20, 90) 
pca.rast.crop <- crop(pca.rast, extent)

# Subset the raster only to include the distribution area (to highlight it)
r2 <- crop(pca.rast.crop, extent(distr.map))
r3 <- mask(r2, distr.map)

# Empty plot to add the graticule of coordinates
empty <- raster(crs=clim.trans)
values(empty) <- 0
empty.crop <- crop(empty, extent)

# Plot PC1 as red, PC2 as green, PC3 as blue
pdf("4-Downstream_Analyses/plots/GDM_fivepcs_RPC1x_GPC2sd_BPC3bio2_cool_RGBmap.pdf", width = 14)
plot(empty.crop, col=c('white'), legend=F)
plotRGB(pca.rast.crop, r=1, g=2, b=3, bgalpha=0, alpha=180, add=T)
plotRGB(r3, r=1, g=2, b=3, bgalpha=0, add=T)
plot(distr.map, lwd=1, add=T)
#points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()
```

## 2. Effects of selection on the genetic structure of populations

To explore this aspect we will first run analyses comparing structure recovered from selected vs neutral loci.

### Extract Neutral Loci



### PCA of candidates

Run PCA with plink on genomics-a
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/
# I manually created a file listing samples to remove for plink with following format:
# c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

plink_1.9 --vcf ll_allsamples_total_intersect_candidate_fivepcs_windows_topsnps.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove samplestoremove.txt \
 --pca --out fivepcs_topsnps_pca
```
Copy eigen files to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/fivepcs_topsnps_pca.eigen\* Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
Prepare R
```{R}
library(tidyverse)
library(RColorBrewer)
library(viridis)
```
plot PCA
```{R}
# import eigen vec and val
pca <- read_table2("4-Downstream_Analyses/tables/fivepcs_topsnps_pca.eigenvec",
                   col_names = FALSE)
eigenval <- scan("4-Downstream_Analyses/tables/fivepcs_topsnps_pca.eigenval")

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

# plot PC1 PC2
pdf(file = paste0("4-Downstream_Analyses/plots/fivepcs_topsnps_pca1_pca2.pdf"),
    width = 8,
    height = 8)

ggplot() + 
  geom_point(data=pca, aes(PC1, PC2), size = 4, pch=21, bg=loc, col="gray32") +
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  scale_color_manual(values=cols) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

dev.off()

# plot PC1 PC3
pdf(file = paste0("4-Downstream_Analyses/plots/fivepcs_topsnps_pca1_pca3.pdf"),
    width = 8,
    height = 8)

ggplot() + 
  geom_point(data=pca, aes(PC1, PC3), size = 4, pch=21, bg=loc, col="gray32") +
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  scale_color_manual(values=cols) +
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

dev.off()
```

### PCA of neutral regions

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

### RDA of candidates

Extract RAW file of topsnps on genomics-a using PLINK
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/
# I manually created a file listing samples to remove for plink with following format:
# c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

plink_1.9 --vcf ll_allsamples_total_intersect_candidate_fivepcs_windows_topsnps.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove samplestoremove.txt --geno 0 \
 --recode A --out fivepcs_topsnps
```
Copy to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/fivepcs_topsnps.raw Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
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
varmatrix <- read_tsv("2-Prepare_Environmental_Data/WorldClim_table_persample.tsv") %>% column_to_rownames(., var="sample")
snowdata <- read_tsv("2-Prepare_Environmental_Data/Snow_table_persample.tsv") %>% column_to_rownames(., var="sample")
env.predictors <- cbind(varmatrix, snowdata)
env.predictors <- env.predictors[, c('bio2', 'bio3', 'snow_days')]

# Load GenoType Data - choose 
gt_data <- read.PLINK("4-Downstream_Analyses/tables/fivepcs_topsnps.raw")
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
# choose file name based on input
pdf(file = paste0("4-Downstream_Analyses/plots/candidate_allsamples_rda1_rda2.pdf"),
    width = 8,
    height = 8)

plot(rda, type="n", scaling=3)
#points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num]) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1)     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

dev.off()

# RDA1 v RDA3
pdf(file = paste0("4-Downstream_Analyses/plots/candidate_allsamples_rda1_rda3.pdf"),
    width = 8,
    height = 8)

plot(rda, type="n", scaling=3, choices=c(1,3))
#points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3, choices=c(1,3))           # the SNPs
points(rda, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num], choices=c(1,3)) # the wolves
text(rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))     # the predictors
#legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg) # legend

dev.off()

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

### PIE CHARTS of each population's allele frequency of inversion

```{R}
allele.freq.table <- data.frame()

for (n in 1:NROW(pca)){
  SAMPLE <- pca$ind[n]
  if (pca$PC1[n] < -0.05){
   INV1 <- 2
   INV2 <- 0
  } else if (pca$PC1[n] > 0.05){
   INV1 <- 0
   INV2 <- 2
  } else if (pca$PC1[n] > -0.05 & pca$PC1[n] < 0.05){
   INV1 <- 1
   INV2 <- 1
  }
  row <- data.frame(sample=SAMPLE, allele.1=INV1, allele.2=INV2)
  allele.freq.table <- rbind(allele.freq.table, row)
}

# Add population
pop <- rep(NA, length(pca$ind))
pop[grep("ca", pca$ind)] <- "Caucasus"
pop[grep("po", pca$ind)] <- "NE-Poland"
pop[grep("ba", pca$ind)] <- "Balkans"
pop[grep("cr", pca$ind)] <- "Carpathians"
pop[grep("no", pca$ind)] <- "Norway"
pop[grep("ka", pca$ind)] <- "Mongolia"
pop[grep("ki", pca$ind)] <- "Kirov"
pop[grep("la", pca$ind)] <- "Latvia"
pop[grep("to", pca$ind)] <- "Mongolia"
pop[grep("og", pca$ind)] <- "Mongolia"
pop[grep("tu", pca$ind)] <- "Tuva"
pop[grep("ur", pca$ind)] <- "Urals"
pop[grep("vl", pca$ind)] <- "Primorsky Krai"
pop[grep("ya", pca$ind)] <- "Yakutia"
allele.freq.table <- cbind(allele.freq.table, pop)

for (n in 1:length(unique(allele.freq.table$pop))){
  POP <- unique(allele.freq.table$pop)[n]
  pop.table <- subset(allele.freq.table, pop == POP)
  
  
  
  al1.freq <- (sum(pop.table$allele.1)) / (2*NROW(pop.table))
  al2.freq <- (sum(pop.table$allele.2)) / (2*NROW(pop.table))
  pie.table <- rbind(data.frame(group="allele.1", value=al1.freq), data.frame(group="allele.2", value=al2.freq))
  
  pie <- ggplot(pie.table, aes(x="", y=value, fill=group)) +
   geom_bar(stat="identity", width=1, color="black") +
   coord_polar("y", start=0) +
   scale_fill_manual(values=c("dodgerblue3", "goldenrod3")) + 
   theme_void() + # remove background, grid, numeric labels
   theme(legend.position="none")
  pie
  ggsave(pie, filename = paste0("4-Downstream_Analyses/plots/inversion_",POP,"_pie.pdf"),  bg = "transparent")

}
```

## 3. Functional annotation of candidate regions

Now we want to check what genes are involved in the adaptive process. To do so, I can intersect the candidate windows with the functional annotation file of my reference genome. The extracted list of genes can be analyzed to discover which functional groups the genes belong to and even test for enrichment of particular functions.

To extract the genes present in the candidate regions:
```{bash}
# On genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/Annotation

# bedtools to intersect bed with gff3 and get a bed of overlapping gene regions
bedtools intersect -wb \
 -a ../Intersect/total_intersect_candidate_windows.bed \
 -b /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.gff3 |
 grep "gbkey=Gene" | cut -f1,2,3,12 \
 > total_intersect_candidate_windows_intersect_genes.bed
# a total of 1018 regions were found

# extract list of genes
cut -d'-' -f2 total_intersect_candidate_windows_intersect_genes.bed | 
 cut -d';' -f1 | sort -u \
 > total_intersect_candidate_windows_intersect_genes_list.txt
# a total of 895 gene IDs were found
# 110 of them have a LOC identifying number from the annotation
# -- CHECK what they mean -- #

# extract list of genes for snow_days and bio2 - TOP ranking variables in GDM
for var in bio2 snow_days
 do
  bedtools intersect -wb \
   -a total_intersect_candidate_windows_intersect_genes.bed \
   -b ../Intersect/${var}_intersect_candidate_windows.bed |
   cut -d'-' -f2 | cut -d';' -f1 | sort -u \
   > ${var}_candidate_windows_intersect_genes_list.txt
done

# extract list of genes in the inversion
bedtools intersect -wb \
 -a ../Intersect/inversion.bed \
 -b /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.gff3 |
 grep "gbkey=Gene" | cut -f1,2,3,12 | cut -d'-' -f2 | cut -d';' -f1 | sort -u \
 > inversion_genes_list.txt
```
Download gene list to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/Annotation/total_intersect_candidate_windows_intersect_genes_list.txt Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/

scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/Annotation/bio2_candidate_windows_intersect_genes_list.txt Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/

scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/Annotation/snow_days_candidate_windows_intersect_genes_list.txt Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/

scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/Annotation/inversion_genes_list.txt Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
Input the list of gene IDs into Panther to make a comparison using the Felis catus annotation for GO-term representation
