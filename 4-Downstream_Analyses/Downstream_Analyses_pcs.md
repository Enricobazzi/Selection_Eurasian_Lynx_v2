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

### GDM

Generalized Dissimilarity Modeling (GDM) will be run to analyze how environmental pressure gradients are affecting genetics differentially in space.

These analyses are based on a mix of these two working examples: <https://github.com/pgugger/LandscapeGenomics/blob/master/2017/Exercise4.md#generalized-dissimilarity-modeling-gdm> <https://www.mountainmanmaier.com/software/pop_genom/#generalized-dissimilarity-modeling>

To get the SNPs in a format readable in R (plink's RAW and then transform):

```{bash}
# on genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/
# get VCF of topsnps only
bedtools intersect -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
 -b Intersect/total_intersect_candidate_fivepcs_windows_topsnps.bed \
 > VCF/ll_allsamples_total_intersect_candidate_fivepcs_windows_topsnps.vcf

cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF

plink_1.9 --vcf ll_allsamples_total_intersect_candidate_fivepcs_windows_topsnps.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove samplestoremove.txt --geno 0 \
 --recode A --out fivepcs_topsnps
```
Download to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/fivepcs_topsnps.raw Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
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
library(factoextra)
library(adegenet)
library(reshape2)
library(vegan)

RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[mapCells] <- residMap
  return(list(max(residMap), rast))
}
```
Import data
```{R}
# Samples to be excluded
samples_to_remove <- c("c_ll_ba_0216","c_ll_ba_0233","c_ll_cr_0211","h_ll_ba_0214","h_ll_ba_0215")

# SNP data
# import RAW of neutral snps
adegen_neutral_raw <-read.PLINK("4-Downstream_Analyses/tables/intergenic_neutral_snps.raw")
neutral_raw <- data.frame(as.matrix(adegen_neutral_raw))

# import RAW of candidates
adegen_candidate_raw <-read.PLINK("4-Downstream_Analyses/tables/fivepcs_topsnps.raw")
candidate_raw <- data.frame(as.matrix(adegen_candidate_raw))

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
# neutral
neutral_snp.dist <- dist(neutral_raw, diag = T, upper = T)   # generate Euclidean distance matrix
neutral_snp.dist.1 <- neutral_snp.dist/(max(neutral_snp.dist))  #rescale by dividing by max value
neutral_snp.dist.1 <- as.matrix(neutral_snp.dist.1)
neutral_snp.dist.1 <- cbind(ID = rownames(neutral_snp.dist.1), neutral_snp.dist.1)  # make the row names an actual column named "ID"
rownames(neutral_snp.dist.1) <- NULL  #remove prior row names

# candidate
candidate_snp.dist <- dist(candidate_raw, diag = T, upper = T)   # generate Euclidean distance matrix
candidate_snp.dist.1 <- candidate_snp.dist/(max(candidate_snp.dist))  #rescale by dividing by max value
candidate_snp.dist.1 <- as.matrix(candidate_snp.dist.1)
candidate_snp.dist.1 <- cbind(ID = rownames(candidate_snp.dist.1), candidate_snp.dist.1)  # make the row names an actual column named "ID"
rownames(candidate_snp.dist.1) <- NULL  #remove prior row names
```
Run GDM. Check relative importance of different variables on shaping genetic dissimilarity between samples. 
```{R}
# Prepare input
neutral_gdm.input <- gdm::formatsitepair(bioData=neutral_snp.dist.1, bioFormat=3, 
                            predData=clim.points, siteColumn="ID", 
                            XColumn="x", YColumn="y")

candidate_gdm.input <- gdm::formatsitepair(bioData=candidate_snp.dist.1, bioFormat=3, 
                            predData=clim.points, siteColumn="ID", 
                            XColumn="x", YColumn="y")

# Run GDM
neutral_gdm <- gdm(neutral_gdm.input, geo = T, splines = NULL, knots = NULL)
candidate_gdm <- gdm(candidate_gdm.input, geo = T, splines = NULL, knots = NULL)

# Print summary
summary(neutral_gdm)
summary(candidate_gdm)

# Overall power of GDM analysis
#   neutral : 55.3799 | candidate : 62.97829 
neutral_gdm$explained
candidate_gdm$explained

# Calculate each variable's importance
neutral_gdm.importance <- gdm.varImp(neutral_gdm.input, geo=T, nPerm=100, parallel=TRUE, cores=8)
candidate_gdm.importance <- gdm.varImp(candidate_gdm.input, geo=T, nPerm=100, parallel=TRUE, cores=8)

# Create data frame for plot - full model importance values
neutral_gdm.imp.df <- data.frame(importance=sort(neutral_gdm.importance[[2]][,1], decreasing=T)) %>% 
  rownames_to_column(., var="variable")
candidate_gdm.imp.df <- data.frame(importance=sort(candidate_gdm.importance[[2]][,1], decreasing=T)) %>% 
  rownames_to_column(., var="variable")

write.table(x = neutral_gdm.imp.df,
             file = paste0("4-Downstream_Analyses/tables/neutral_gdm_importance.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")
write.table(x = candidate_gdm.imp.df,
             file = paste0("4-Downstream_Analyses/tables/candidate_gdm_importance.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

# Plot importance in a barplot ordered from most to least important
p <- ggplot(data=neutral_gdm.imp.df, aes(x=variable, y=importance)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=c(rev(neutral_gdm.imp.df$variable))) +
  theme_minimal() +
  coord_flip()

pdf(file = paste0("4-Downstream_Analyses/plots/GDM_neutral_importance.pdf"),
    width = 8,
    height = 8)
p
dev.off()

p <- ggplot(data=candidate_gdm.imp.df, aes(x=variable, y=importance)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=c(rev(candidate_gdm.imp.df$variable))) +
  theme_minimal() +
  coord_flip()

pdf(file = paste0("4-Downstream_Analyses/plots/GDM_fivepcs_candidate_importance.pdf"),
    width = 8,
    height = 8)
p
dev.off()

# Plot the effect function of each variable, with uncertainty in estimation (confidence interval?)
pdf(file = paste0("4-Downstream_Analyses/plots/GDM_fivepcs_neutral_effect_uncertainty.pdf"),
    width = 8,
    height = 8)
plotUncertainty(neutral_gdm.input, sampleSites=0.70, bsIters=100, geo=T, plot.layout=c(2,4))
dev.off()

pdf(file = paste0("4-Downstream_Analyses/plots/GDM_fivepcs_candidate_effect_uncertainty.pdf"),
    width = 8,
    height = 8)
plotUncertainty(candidate_gdm.input, sampleSites=0.70, bsIters=100, geo=T, plot.layout=c(2,4))
dev.off()

neutral_gdm.imp.df$variable[1]

neutral_gdm.imp.df$importance[1]/sum(neutral_gdm.imp.df$importance)*100

# Plot heat(?) of relative importance for each variable - USE INSTEAD OF GDM.IMP??
# meant as maximum value of Ispline
plot(neutral_gdm)
plot(candidate_gdm)
# extract isplines from GDMs
neutral_Splines <- isplineExtract(neutral_gdm) # extract spline data for custom plotting
candidate_Splines <- isplineExtract(candidate_gdm) # extract spline data for custom plotting

neutral_importance_base <- data.frame(neutral_Splines$x)
neutral_importance_curve <- data.frame(neutral_Splines$y)
candidate_importance_base <- data.frame(candidate_Splines$x)
candidate_importance_curve <- data.frame(candidate_Splines$y)

#PLOT???#

# create maximum spline (importance) matrix
neutral_importance <- neutral_importance_curve %>% summarise_if(is.numeric, max)
candidate_importance <- candidate_importance_curve %>% summarise_if(is.numeric, max)
neutral_importance_df <- data.frame(importance=t(neutral_importance)) %>% rownames_to_column("variable")
neutral_importance_df <- neutral_importance_df[order(-neutral_importance_df$importance),]

candidate_importance_df <- data.frame(importance=t(candidate_importance)) %>% rownames_to_column("variable")
candidate_importance_df <- candidate_importance_df[order(-candidate_importance_df$importance),]

p <- ggplot(data=filter(neutral_importance_df, importance > 0.01), aes(x=variable, y=importance)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=c(rev(filter(neutral_importance_df, importance > 0.01)[,1]))) +
  #scale_y_continuous(limits = c(0,1))+
  theme_minimal(base_size = 24) +
  theme(axis.title.y = element_text(size = 0) )+
  coord_flip()
pdf(file = paste0("4-Downstream_Analyses/plots/GDM_neutral_isplineimportance.pdf"),
    width = 6,
    height = 8)
p
dev.off()

p <- ggplot(data=filter(candidate_importance_df, importance > 0.01), aes(x=variable, y=importance)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=c(rev(filter(candidate_importance_df, importance > 0.01)[,1]))) +
  #scale_y_continuous(limits = c(0,1))+
  theme_minimal(base_size = 24) +
  theme(axis.title.y = element_text(size = 0) )+
  coord_flip()
pdf(file = paste0("4-Downstream_Analyses/plots/GDM_candidate_isplineimportance.pdf"),
    width = 6,
    height = 8)
p
dev.off()


importance <- rbind(neutral_importance,candidate_importance)
row.names(importance)[c(1,2)] <- c("neutral", "candidate")
importance <- importance %>% rownames_to_column(var="set")
importance.1 <- melt(importance) %>% filter(value>0.1)
importance.1[is.na(importance.1)] = 0
# plot importance matrix
pdf(file = paste0("4-Downstream_Analyses/plots/GDM_isplineimportance.pdf"),
    width = 8,
    height = 8)
ggplot(importance.1, aes(set, variable, fill= value)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_minimal()
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
neutral_clim.trans <- gdm.transform(neutral_gdm, clim.layer)
candidate_clim.trans <- gdm.transform(candidate_gdm, clim.layer)

# Write rasters
writeRaster(neutral_clim.trans, "4-Downstream_Analyses/tables/neutral_gdm_raster.tif", 
            format="GTiff", overwrite=TRUE)
writeRaster(candidate_clim.trans, "4-Downstream_Analyses/tables/candidate_gdm_raster.tif", 
            format="GTiff", overwrite=TRUE)
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
Plot RGB map of 2 most affecting variables (bio2 and snow_days)
```{R}
extent <- c(-10, 180, 20, 80)

empty <- raster(crs=neutral_clim.trans)
values(empty) <- 0
empty.crop <- crop(empty, extent)

# check which layers are which variable:
candidate_clim.trans[[4]]@data@names # bio2
# crop map of bio2
bio2_map.crop <- crop(candidate_clim.trans[[4]], extent)
bio2_map_r2 <- crop(bio2_map.crop, extent(distr.map))
bio2_map_r3 <- mask(bio2_map_r2, distr.map)

## HAVE TO MODIFY THIS IN ILLUSTRATOR OR EQUIVALENT AS DISTRIBUTION MAP ENDS UP SHIFTED I DON'T KNOW WHYY!?!?!11!? ##
pdf("4-Downstream_Analyses/plots/GDM_candidate_bio2_map.pdf", width = 12.5)
plot(empty.crop, col=c('white'), legend=F)
plot(bio2_map.crop, alpha=0.7, col=rev(brewer.pal(8, "YlOrRd")), legend=F, add=T)
plot(bio2_map_r3, col=rev(brewer.pal(8, "YlOrRd")), add=T)
plot(distr.map, lwd=1, add=T)
points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()

# check which layers are which variable:
candidate_clim.trans[[9]]@data@names # snow_days
# crop map of bio2
snow_days_map.crop <- crop(candidate_clim.trans[[9]], extent)
snow_days_map_r2 <- crop(snow_days_map.crop, extent(distr.map))
snow_days_map_r3 <- mask(snow_days_map_r2, distr.map)

## HAVE TO MODIFY THIS IN ILLUSTRATOR OR EQUIVALENT AS DISTRIBUTION MAP ENDS UP SHIFTED I DON'T KNOW WHYY!?!?!11!? ##
pdf("4-Downstream_Analyses/plots/GDM_candidate_snow_days_map.pdf", width = 12.5)
plot(empty.crop, col=c('white'), legend=F)
plot(snow_days_map.crop, alpha=0.7, col=rev(brewer.pal(8, "YlOrRd")), legend=F, add=T)
plot(snow_days_map_r3, col=rev(brewer.pal(8, "YlOrRd")), add=T)
plot(distr.map, lwd=1, add=T)
points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()


############################################################################################################
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
neutral_clim.trans <- gdm.transform(neutral_gdm, clim.layer)
candidate_clim.trans <- gdm.transform(candidate_gdm, clim.layer)

# Extract values from the transformed raster
neutral_clim.rast <- na.omit(getValues(neutral_clim.trans))
candidate_clim.rast <- na.omit(getValues(candidate_clim.trans))

# run PCA on raster values
neutral_pca <- prcomp(neutral_clim.rast)
neutral_pca_table <- data.frame(neutral_pca$rotation) %>% mutate_if(is.numeric, round, digits=3) %>% 
  rownames_to_column(var="predictor")
write.table(x = neutral_pca_table,
             file = paste0("4-Downstream_Analyses/tables/gdm_neutral_pca_table.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

#                  PC1           PC2          PC3          PC4           PC5
# xCoord  9.986830e-01  0.0229120809  0.031788107  0.026929314 -0.0187554932 
# yCoord -2.420250e-02  0.6371006937 -0.189084136  0.321600054 -0.3969175873 
# bio1    1.815996e-03 -0.2252073373  0.169368706  0.016779698  0.1754847956 
# bio2    5.342185e-03 -0.0279138020  0.016938120 -0.018936839  0.0660655704 
# bio3   -1.333081e-02 -0.6290162208  0.165526783  0.434945367 -0.5962496602 
# bio5    1.042126e-02 -0.1332319052  0.125607230 -0.119025330  0.3145769849 
# bio10   5.548369e-03 -0.1636302589  0.090230888 -0.016916564  0.2851002156 
# bio12  -2.126604e-05  0.0008019347  0.001898323 -0.004635841 -0.0036819107 
# bio13   4.066864e-04  0.0071598856  0.021282124 -0.128173733 -0.0949470053 
# bio15   4.034265e-02 -0.3186886504 -0.929196824 -0.125720394 -0.0009901809 
# bio18   8.353938e-03 -0.0166986931  0.138140770 -0.811986311 -0.5120167447

candidate_pca <- prcomp(candidate_clim.rast)
candidate_pca_table <- data.frame(candidate_pca$rotation) %>% mutate_if(is.numeric, round, digits=3) %>% 
  rownames_to_column(var="predictor")
write.table(x = candidate_pca_table,
             file = paste0("4-Downstream_Analyses/tables/gdm_candidate_pca_table.tsv"),
             quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

#                         PC1         PC2          PC3         PC4         PC5 
# xCoord          0.994665893  0.05137837  0.082704592 -0.00741857 -0.01168617 
# yCoord         -0.027428332  0.36517133  0.071363979  0.50047068 -0.46300817 
# bio1            0.002101993 -0.10979327  0.035121077 -0.02295548  0.04938283 
# bio2            0.093201157 -0.30745768 -0.921091013  0.19020439 -0.05558361 
# bio3           -0.009612567 -0.32107801  0.039958275 -0.63082480 -0.35286185 
# bio15           0.026526467 -0.09031071 -0.004241555 -0.31700978 -0.17302210 
# bio18           0.009849161 -0.01713787  0.016389079 -0.04504736  0.76243233 
# bio19          -0.007645500  0.01918360 -0.002308274  0.10538952  0.20480199 
# mean_snow_days -0.015654397  0.80344045 -0.369508327 -0.44855106  0.04757753

# Create a raster with PC values
candidate_pca.rast <- predict(candidate_clim.trans, candidate_pca, index=1:5)
neutral_pca.rast <- predict(neutral_clim.trans, neutral_pca, index=1:5)

# Transform the values in the wanted layers to a scale of 255 for RGB colors
neutral_pca.rast[[1]] <- (neutral_pca.rast[[1]]-neutral_pca.rast[[1]]@data@min) / (neutral_pca.rast[[1]]@data@max-neutral_pca.rast[[1]]@data@min)*255
neutral_pca.rast[[2]] <- (neutral_pca.rast[[2]]-neutral_pca.rast[[2]]@data@min) / (neutral_pca.rast[[2]]@data@max-neutral_pca.rast[[2]]@data@min)*255
neutral_pca.rast[[3]] <- (neutral_pca.rast[[3]]-neutral_pca.rast[[3]]@data@min) / (neutral_pca.rast[[3]]@data@max-neutral_pca.rast[[3]]@data@min)*255
neutral_pca.rast[[4]] <- (neutral_pca.rast[[4]]-neutral_pca.rast[[4]]@data@min) / (neutral_pca.rast[[4]]@data@max-neutral_pca.rast[[4]]@data@min)*255
neutral_pca.rast[[5]] <- (neutral_pca.rast[[5]]-neutral_pca.rast[[5]]@data@min) / (neutral_pca.rast[[5]]@data@max-neutral_pca.rast[[5]]@data@min)*255
candidate_pca.rast[[1]] <- (candidate_pca.rast[[1]]-candidate_pca.rast[[1]]@data@min) / (candidate_pca.rast[[1]]@data@max-candidate_pca.rast[[1]]@data@min)*255
candidate_pca.rast[[2]] <- (candidate_pca.rast[[2]]-candidate_pca.rast[[2]]@data@min) / (candidate_pca.rast[[2]]@data@max-candidate_pca.rast[[2]]@data@min)*255
candidate_pca.rast[[3]] <- (candidate_pca.rast[[3]]-candidate_pca.rast[[3]]@data@min) / (candidate_pca.rast[[3]]@data@max-candidate_pca.rast[[3]]@data@min)*255
candidate_pca.rast[[4]] <- (candidate_pca.rast[[4]]-candidate_pca.rast[[4]]@data@min) / (candidate_pca.rast[[4]]@data@max-candidate_pca.rast[[4]]@data@min)*255
candidate_pca.rast[[5]] <- (candidate_pca.rast[[5]]-candidate_pca.rast[[5]]@data@min) / (candidate_pca.rast[[5]]@data@max-candidate_pca.rast[[5]]@data@min)*255

# I was trying to make a cool legend - this crap does not work:
# legend.colors <- na.omit(getValues(pca.rast[[1:3]]))
# colors <- c()
# for (r in 1:NROW(legend.colors)){
#   print(r)
#   red <- legend.colors[r,1]/255
#   green <- legend.colors[r,2]/255
#   blue <- legend.colors[r,3]/255
#   colo <- rgb(red,green,blue,alpha=1)
#   colors <- c(colors,colo)
# }
# pdf("~/Desktop/legend.pdf", width = 12.5)
# fviz_pca_biplot(pca, repel = TRUE, label = "var",
#                 col.var = "black", # Variables color
#                 col.ind = colors  # Individuals color
#                 )
# dev.off()
##############
# unsupervised classification:
# library("randomForest")
# # kmeans:
# v <- getValues(r3[[2:4]])
# i <- which(!is.na(v))
# v <- na.omit(v)
# E <- kmeans(v, 4, iter.max = 100, nstart = 10)
# kmeans_raster <- raster(r3[[2:4]])
# 
# kmeans_raster[i] <- E$cluster
# plot(kmeans_raster, col=rainbow(4))
# 
# ## unsupervised randomForest classification using kmeans
# vx<-v[sample(nrow(v), 10000),]
# rf = randomForest(vx)
# rf_prox <- randomForest(vx,ntree = 1000, proximity = TRUE)$proximity
# 
# E_rf <- kmeans(rf_prox, 4, iter.max = 100, nstart = 10)
# rf <- randomForest(vx,as.factor(E_rf$cluster),ntree = 500)
# rf_raster<- predict(r3[[1:3]],rf)
# plot(rf_raster, col=rainbow(4))


##############

# set the extent of map to plot only interesting part
extent <- c(-10, 180, 20, 80)
neutral_pca.rast.crop <- crop(neutral_pca.rast, extent)
candidate_pca.rast.crop <- crop(candidate_pca.rast, extent)

# Subset the raster only to include the distribution area (to highlight it)
neutral_r2 <- crop(neutral_pca.rast.crop, extent(distr.map))
neutral_r3 <- mask(neutral_r2, distr.map)
candidate_r2 <- crop(candidate_pca.rast.crop, extent(distr.map))
candidate_r3 <- mask(candidate_r2, distr.map)

# Empty plot to add the graticule of coordinates
empty <- raster(crs=neutral_clim.trans)
values(empty) <- 0
empty.crop <- crop(empty, extent)

# Plot ONLY PC1, PC2 and PC3 (~xCoord, ~snow and ~bio2)
candidate_pc2.rast.crop <- candidate_pca.rast.crop
values(candidate_pc2.rast.crop[[5]]) <- 0
candidate_pc2_r2 <- crop(candidate_pc2.rast.crop, extent(distr.map))
candidate_pc2_r3 <- mask(candidate_pc2_r2, distr.map)
# PC1
plot(empty.crop, col=c('white'), legend=F)
plotRGB(candidate_pc2.rast.crop, r=1, g=5, b=5, bgalpha=0, alpha=180, add=T)
plotRGB(candidate_pc2_r3, r=1, g=5, b=5, bgalpha=0, add=T)

# PC2
plot(empty.crop, col=c('white'), legend=F)
plotRGB(candidate_pc2.rast.crop, r=5, g=2, b=5, bgalpha=0, alpha=180, add=T)
plotRGB(candidate_pc2_r3, r=5, g=2, b=5, bgalpha=0, add=T)

pdf("4-Downstream_Analyses/plots/GDM_candidate_PC2_map.pdf", width = 12.5)
plot(empty.crop, col=c('white'), legend=F)
plot(candidate_pc2.rast.crop[[2]], alpha=0.7, col=rev(brewer.pal(8, "YlOrRd")), legend=F, add=T)
plot(candidate_pc2_r3[[2]], col=rev(brewer.pal(8, "YlOrRd")), add=T)
plot(distr.map, lwd=1, add=T)
points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()

# PC3
plot(empty.crop, col=c('white'), legend=F)
plotRGB(candidate_pc2.rast.crop, r=5, g=5, b=3, bgalpha=0, alpha=180, add=T)
plotRGB(candidate_pc2_r3, r=5, g=5, b=3, bgalpha=0, add=T)

pdf("4-Downstream_Analyses/plots/GDM_candidate_PC3_map.pdf", width = 12.5)
plot(empty.crop, col=c('white'), legend=F)
plot(candidate_pc2.rast.crop[[3]], alpha=0.7, col=rev(brewer.pal(4, "YlOrRd")), legend=F, add=T)
plot(candidate_pc2_r3[[3]], col=rev(brewer.pal(4, "YlOrRd")), add=T)
plot(distr.map, lwd=1, add=T)
points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()


# Plot PC1 as red, PC2 as green, PC3 as blue
pdf("4-Downstream_Analyses/plots/GDM_neutral_RPC1x_GPC2sd_BPC3bio2_cool_RGBmap.pdf", width = 12.5)
plot(empty.crop, col=c('white'), legend=F)
plotRGB(neutral_pca.rast.crop, r=1, g=2, b=3, bgalpha=0, alpha=180, add=T)
plotRGB(neutral_r3, r=1, g=2, b=3, bgalpha=0, add=T)
plot(distr.map, lwd=1, add=T)
points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()

pdf("4-Downstream_Analyses/plots/GDM_fivepcs_candidate_RPC1x_GPC2sd_BPC3bio2_cool_RGBmap.pdf", width = 12.5)
plot(empty.crop, col=c('white'), legend=F)
plotRGB(candidate_pca.rast.crop, r=1, g=2, b=3, bgalpha=0, alpha=180, add=T)
plotRGB(candidate_r3, r=1, g=2, b=3, bgalpha=0, add=T)
plot(distr.map, lwd=1, add=T)
points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()

# difference
RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[c(mapCells)] <- residMap
  return(list(max(residMap), rast))
}
env_trns <- data.frame(getValues(neutral_clim.trans)) %>% rownames_to_column("cells") %>% drop_na()
MASK <- neutral_clim.trans$xCoord
MASK[MASK > -999 ] <- NA
plot(MASK)
diff_map <- RGBdiffMap(neutral_clim.rast,candidate_clim.rast,MASK,as.numeric(env_trns$cell))
diff_map.crop <- crop(diff_map[[2]], extent)
diff_map_r2 <- crop(diff_map.crop, extent(distr.map))
diff_map_r3 <- mask(diff_map_r2, distr.map)

# Subset the raster only to include the distribution area (to highlight it)
## HAVE TO MODIFY THIS IN ILLUSTRATOR OR EQUIVALENT AS DISTRIBUTION MAP ENDS UP SHIFTED I DON'T KNOW WHYY!?!?!11!? ##
pdf("4-Downstream_Analyses/plots/GDM_neutral_candidate_diffProcrust_map.pdf", width = 12.5)
plot(empty.crop, col=c('white'), legend=F)
plot(diff_map.crop, alpha=0.7, col=rev(brewer.pal(6, "RdYlBu")), legend=F, add=T)
plot(diff_map_r3, col=rev(brewer.pal(6, "RdYlBu")), add=T)
plot(distr.map, lwd=1, add=T)
points(clim.points$x,clim.points$y, pch=21, lwd=0.7, cex=0.8, col="black", bg=bg[num]) # the lynxes
dev.off()
```

## 2. Effects of selection on the genetic structure of populations

To explore this aspect we will first run analyses comparing structure recovered from selected vs neutral loci.

### Extract Neutral Loci

To extract only neutral SNPs I will filter the VCF file to remove any candidate regions and all of the genes genes.

Genes regions were identified by Dani as follows:
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
```
Remove those regions from the VCF
```{bash}
# on genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF

bedtools subtract -header -a ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
 -b /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.genes.plus5000.bed \
 > ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.vcf
```
From 4983054 whole genome SNPs to 2476154 intergenic SNPs.

Now remove also any candidate of selection - the SNPs from the RDA and the GenWin windows from BayPass
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/
# create a new VCF to substract from sequentially
cp VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.vcf \
 VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf
 
# Subtract GenWin windows sequentially from new VCF
for var in PC1 PC2 PC3 PC4 PC5
 do
  echo "filtering ${var} regions"
  bedtools subtract -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
   -b GenWin/${var}_GenWin_windows_outliers.bed \
   > tmp && mv tmp VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf
done

# Subtract RDA SNPs from new VCF
bedtools subtract -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
 -b RDA/rda_candidate_fivepcs_snps.bed \
 > tmp && mv tmp VCF/ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf
```
We have a total of 2266071 non-selected intergenic SNPs

Now to create a RAW file with plink
```{bash}
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/
# I manually created a file listing samples to remove for plink with following format:
# c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

# prune snps based on ld (VIF < 2) - no missing data - maf 0.05
plink_1.9 --vcf ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# --geno 0 --maf 0.05 \
--remove samplestoremove.txt --indep 100 10 2

# Extract RAW file
plink_1.9 --vcf ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--remove samplestoremove.txt --geno 0 --extract plink.prune.in \
--recode A --out intergenic_neutral_snps
```
Copy to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/intergenic_neutral_snps.raw Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```

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

# plot PC2 PC3
pdf(file = paste0("4-Downstream_Analyses/plots/fivepcs_topsnps_pca2_pca3.pdf"),
    width = 8,
    height = 8)

ggplot() + 
  geom_point(data=pca, aes(PC2, PC3), size = 4, pch=21, bg=loc, col="gray32") +
  coord_equal() + theme_light() + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
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
 --remove samplestoremove_pca.txt --extract <(shuf --random-source=<(yes 42) plink.prune.in | head -n 10000) \
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
library(viridis)
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

# plot PC1 vs PC2
pdf(file = paste0("4-Downstream_Analyses/plots/neutral_pca1_pca2.pdf"),
    width = 8,
    height = 8)

ggplot() + 
  geom_point(data=pca, aes(PC1, PC2), size = 4, pch=21, bg=loc, col="gray32") +
  coord_equal(0.5) + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  
  scale_color_manual(values=cols) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

dev.off()

# plot PC1 vs PC3
pdf(file = paste0("4-Downstream_Analyses/plots/neutral_pca1_pca3.pdf"),
    width = 8,
    height = 8)

ggplot() + 
  geom_point(data=pca, aes(PC1, PC3), size = 4, pch=21, bg=loc, col="gray32") +
  coord_equal(0.5) + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  scale_color_manual(values=cols) +
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

dev.off()
```

### Discriminant analysis of Principal Components (DAPC)

Generate RAW files with plink
```{bash}
# on genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/

# neutral raw data for ade4 dudi.pca
plink_1.9 --vcf ll_wholegenome_LyCa_ref.sorted.filter7.intergenic.noselection.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove samplestoremove_pca.txt --extract <(shuf --random-source=<(yes 42) plink.prune.in | head -n 1500) \
 --recode A --out intergenic_neutral_snps_dudi

# candidate raw data for ade4 dudi.pca
plink_1.9 --vcf ll_allsamples_total_intersect_candidate_fivepcs_windows_topsnps.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove samplestoremove.txt --geno 0 \
 --recode A --out fivepcs_topsnps_dudi
```
Download to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/VCF/\*dudi.raw Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```

```{R}
# prepare R
library(tidyverse)
library(adegenet)
library(vegan)
library(RColorBrewer)
library(viridis)
library(factoextra)

# Load GenoType Data - choose 
# import RAW of neutral snps
adegen_neutral_raw <-read.PLINK("4-Downstream_Analyses/tables/intergenic_neutral_snps.raw")
neutral_raw <- data.frame(as.matrix(adegen_neutral_raw))

# import RAW of candidates
adegen_candidate_raw <-read.PLINK("4-Downstream_Analyses/tables/fivepcs_topsnps.raw")
candidate_raw <- data.frame(as.matrix(adegen_candidate_raw))

# check data using pca
pca_neutr <- dudi.pca(neutral_raw, scale = TRUE, scan = FALSE, nf = 3)
pca_cand <- dudi.pca(candidate_raw, scale = TRUE, scan = FALSE, nf = 3)
fviz_pca_ind(pca_neutr, axes = c(1,2),
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
             #repel = TRUE     # Avoid text overlapping
             )
fviz_pca_ind(pca_cand, axes = c(1,2),
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
             #repel = TRUE     # Avoid text overlapping
             )

# DAPC of neutral
grp <- find.clusters(adegen_neutral_raw, max.n.clust=20)
grps <- grp$grp
neutral_grps <- grp$grp

xval <- xvalDapc(adegen_neutral_raw, grps, n.pca.max = 200,
                 result = "overall",
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)

DAPC <- dapc(adegen_neutral_raw, grps)

DAPC$posterior
myCol <- c(brewer.pal(12,"Paired")[5],"#A035AF","#B8860b",viridis_pal()(5)[1],brewer.pal(12,"Paired")[3])
scatter(DAPC, col=myCol, scree.da=FALSE, xax = 1, yax = 4,
        cell=1.5, cex=2, bg="white",cstar=1)

scatter(DAPC)

# DAPC of candidates
grp <- find.clusters(adegen_candidate_raw, n.pca = 103, max.n.clust = 20)
kstat <- grp$Kstat

grps <- grp$grp
candidate_grps <- grp$grp

xval <- xvalDapc(adegen_candidate_raw, candidate_grps, n.pca.max = 200,
                 result = "overall",
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)

DAPC <- dapc(adegen_candidate_raw, candidate_grps, n.pca = xval$DAPC$n.pca, n.da = xval$DAPC$n.da)
DAPC$posterior
scatter(DAPC)
myCol <- c(brewer.pal(n = 9, name = "Set3"))
scatter(DAPC, col=myCol, scree.da=FALSE, xax = 1, yax = 3,
        cell=1.5, cex=2, bg="white",cstar=1)

set.seed(4)
contrib <- loadingplot(DAPC$var.contr, axis=1,
                       thres=.008, lab.jitter=1)


DAPC$var.contr

igraph::compare(candidate_grps,neutral_grps, method="adjusted.rand")

```


### Co-Inertia of Neutral vs Candidate loci

Run coinertia
```{R}
library(tidyverse)
library(ade4)
library(adegenet)
library(made4)
library(RColorBrewer)
library(viridis)

# import RAW of neutral snps
adegen_neutral_raw <-read.PLINK("4-Downstream_Analyses/tables/intergenic_neutral_snps.raw")
neutral_raw <- data.frame(as.matrix(adegen_neutral_raw))

# import RAW of candidates
adegen_candidate_raw <-read.PLINK("4-Downstream_Analyses/tables/fivepcs_topsnps.raw")
candidate_raw <- data.frame(as.matrix(adegen_candidate_raw))

# remove samples not in common
samples_to_remove <- c(setdiff(row.names(candidate_raw), row.names(neutral_raw)))
candidate_raw <- candidate_raw[!(row.names(candidate_raw) %in% samples_to_remove),]

# using made4

# coinertia:
cia1 <- cia(t(neutral_raw), t(candidate_raw), cia.nf=3, cia.scan=FALSE, nsc=TRUE)
# calculate RV = between 0 (completely different) and 1 (the same)
RV=cia1$coinertia$RV
# using 3 dimensions (plot in 3D -> still need to work on it)
cia_df <- data.frame(sample=(rownames(cia1$coinertia$mX)), 
                       start_x=(cia1$coinertia$mX$NorS1), start_y=(cia1$coinertia$mX$NorS2), start_z=(cia1$coinertia$mX$NorS3),
                       end_x=(cia1$coinertia$mY$NorS1), end_y=(cia1$coinertia$mY$NorS2), end_z=(cia1$coinertia$mY$NorS3))
# using 2 dimensions
cia_df <- data.frame(sample=(rownames(cia1$coinertia$mX)), 
                       start_x=(cia1$coinertia$mX$NorS1), start_y=(cia1$coinertia$mX$NorS2),
                       end_x=(cia1$coinertia$mY$NorS1), end_y=(cia1$coinertia$mY$NorS2))

loc <- rep(NA, length(cia_df$sample))
loc[grep("ca", cia_df$sample)] <- "#B8860b"
loc[grep("po", cia_df$sample)] <- viridis_pal()(5)[3]
loc[grep("ba", cia_df$sample)] <- "#A035AF"
loc[grep("cr", cia_df$sample)] <- brewer.pal(12,"Paired")[9]
loc[grep("no", cia_df$sample)] <- viridis_pal()(5)[2]
loc[grep("ka", cia_df$sample)] <- brewer.pal(12,"Paired")[7]
loc[grep("ki", cia_df$sample)] <- viridis_pal()(5)[1]
loc[grep("la", cia_df$sample)] <- brewer.pal(12,"Paired")[3]
loc[grep("to", cia_df$sample)] <- brewer.pal(12,"Paired")[7]
loc[grep("og", cia_df$sample)] <- brewer.pal(12,"Paired")[7]
loc[grep("tu", cia_df$sample)] <- brewer.pal(12,"Paired")[8]
loc[grep("ur", cia_df$sample)] <- "#0F4909"
loc[grep("vl", cia_df$sample)] <- brewer.pal(12,"Paired")[5]
loc[grep("ya", cia_df$sample)] <- brewer.pal(12,"Paired")[6]
cia_df <- data.frame(cia_df, color=(loc))

# distance in 3D
cia_df$distance <- apply(cia_df[,c(2:7)], 1, function(x) dist(matrix(x, nrow = 1, byrow = TRUE)))
# distance in 2D
cia_df$distance <- apply(cia_df[,c(2:5)], 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
# standardized distance
#cia_df$distance <- abs(cia_df$distance - mean(cia_df$distance))/(sd(cia_df$distance))

distances <- c()
for (i in 1:nrow(cia_df)){
  distance <- sqrt(((cia_df[i,4]-cia_df[i,2])**2) + ((cia_df[i,5]-cia_df[i,3])**2))
  distances[i] <- distance
}
cia_df$distance2 <- distances
#cia_df$distance2 <- abs((cia_df$distance2 - mean(cia_df$distance2))/(sd(cia_df$distance2)))

# plot coinertia calculated trajectories, with point size proportional to standardized distance between neutral and candidate
# position in coinertia space
pdf(file = paste0("4-Downstream_Analyses/plots/coinertia_axis1_axis2_stdist.pdf"),
    width = 16,
    height = 6)
ggplot() + 
  geom_point(data=cia_df, mapping=aes(x=start_x, y=start_y), size=0.5, fill=loc) +
  geom_segment(data=cia_df, mapping=aes(x=start_x, y=start_y, xend=end_x, yend=end_y), 
               arrow = arrow(length = unit(0.25,"cm")), size=0.4, color="darkgrey", alpha=0.8) + 
  geom_point(data=cia_df, mapping=aes(x=end_x, y=end_y), size=distances*80, shape=21, fill=loc, alpha=0.6) +
  scale_y_continuous(limits=c(-1.25,0.65))+
theme_minimal()
dev.off()

# add population average shift arrows
pops <- c("ba", "ca", "cr", "ki", "la", "no", "po", "tu", "ur", "vl", "ya")
avg_arrows_df <- data.frame()
for (n in 1:length(pops)){
 pop=pops[n]
 rows <- grep(pop, cia_df$sample)
 pop_cia <- cia_df[c(rows),]
 pop_arrow <- data.frame(popu=pop, start_x=mean(pop_cia$start_x),start_y=mean(pop_cia$start_y),
                         end_x=mean(pop_cia$end_x),end_y=mean(pop_cia$end_y))
 distance <- sqrt(((pop_arrow[1,4]-pop_arrow[1,2])**2) + ((pop_arrow[1,5]-pop_arrow[1,3])**2))
 pop_arrow$distance2 <- distance
 avg_arrows_df <- data.frame(rbind(avg_arrows_df, pop_arrow))
}

ggplot() + 
  geom_segment(data=avg_arrows_df, mapping=aes(x=start_x, y=start_y, xend=end_x, yend=end_y), 
               arrow = arrow(), size=2, color=(unique(cia_df$color)[c(1:3,5:12)]), alpha=0.8) + 
  #scale_y_continuous(limits=c(-1.25,0.65))+
theme_minimal()

```
Run coinertia on allele frequencies
```{R}
# CALCULATE NEUTRAL ALLELE FREQUENCIES
# import RAW of neutral snps
adegen_neutral_raw <-read.PLINK("4-Downstream_Analyses/tables/intergenic_neutral_snps.raw")
neutral_raw <- data.frame(as.matrix(adegen_neutral_raw))
# calculate neutral allele frequencies - no mongolia
pops <- c("ba", "ca", "cr", "ki", "la", "no", "po", "tu", "ur", "vl", "ya")
af_dataset <- data.frame()

for (n in 1:length(pops)){
 pop=pops[n]

 rows <- grep(pop, ROWNAMES(neutral_raw))

 pop_alleles <- neutral_raw[c(rows),]

 pop_afs <- c()
 for (i in 1:ncol(pop_alleles)){
  locus <- pop_alleles[,i]
  tot_alleles <- nrow(pop_alleles)*2
  locus_f <- sum(locus) / tot_alleles
  pop_afs <- c(pop_afs, locus_f)
  }

af_dataset <- rbind(af_dataset, pop_afs)
row.names(af_dataset)[n] <- pop  
}

# add mongolia allele frequencies
mong <- c("ka", "og", "to")
rows <- c()
for (p in mong){
  app <- grep(p,ROWNAMES(neutral_raw))
  rows <- c(rows, app)
}
pop_alleles <- neutral_raw[c(rows),]
pop_afs <- c()
for (i in 1:ncol(pop_alleles)){
  locus <- pop_alleles[,i]
  tot_alleles <- nrow(pop_alleles)*2
  locus_f <- sum(locus) / tot_alleles
  pop_afs <- c(pop_afs, locus_f)
}
pop="mo"
n=12
af_dataset <- rbind(af_dataset, pop_afs)
row.names(af_dataset)[n] <- pop  
colnames(af_dataset) <- c(colnames(neutral_raw))

# CALCULATE CANDIDATE ALLELE FREQUENCIES
# import RAW of candidates
adegen_candidate_raw <-read.PLINK("4-Downstream_Analyses/tables/fivepcs_topsnps.raw")
candidate_raw <- data.frame(as.matrix(adegen_candidate_raw))
# calculate candidate allele frequencies - no mongolia
pops <- c("ba", "ca", "cr", "ki", "la", "no", "po", "tu", "ur", "vl", "ya")
af_candidate_dataset <- data.frame()

for (n in 1:length(pops)){
 pop=pops[n]

 rows <- grep(pop, ROWNAMES(candidate_raw))

 pop_alleles <- candidate_raw[c(rows),]

 pop_afs <- c()
 for (i in 1:ncol(pop_alleles)){
  locus <- pop_alleles[,i]
  tot_alleles <- nrow(pop_alleles)*2
  locus_f <- sum(locus) / tot_alleles
  pop_afs <- c(pop_afs, locus_f)
  }

af_candidate_dataset <- rbind(af_candidate_dataset, pop_afs)
row.names(af_candidate_dataset)[n] <- pop  
}

# add mongolia allele frequencies
mong <- c("ka", "og", "to")
rows <- c()
for (p in mong){
  app <- grep(p,ROWNAMES(candidate_raw))
  rows <- c(rows, app)
}
pop_alleles <- candidate_raw[c(rows),]
pop_afs <- c()
for (i in 1:ncol(pop_alleles)){
  locus <- pop_alleles[,i]
  tot_alleles <- nrow(pop_alleles)*2
  locus_f <- sum(locus) / tot_alleles
  pop_afs <- c(pop_afs, locus_f)
}
pop="mo"
n=12
af_candidate_dataset <- rbind(af_candidate_dataset, pop_afs)
row.names(af_candidate_dataset)[n] <- pop  
colnames(af_candidate_dataset) <- c(colnames(candidate_raw))

# RUN COINERTIA
# Run dudi.pca on 2 sets
library("factoextra")
pca_neutr <- dudi.pca(af_dataset, scale = TRUE, scan = FALSE, nf = 3)
pca_cand <- dudi.pca(af_candidate_dataset, scale = TRUE, scan = FALSE, nf = 3)

fviz_pca_ind(pca_cand,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
fviz_pca_ind(pca_neutr, axes = c(1,2),
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )


coin1 <- coinertia(pca_neutr, pca_cand, scan = FALSE, nf = 2)

summary(coin1)

pdf(file = paste0("~/Desktop/coin1.pdf"),
    width = 8,
    height = 8)
plot(coin1, clab=0, nlab=0)
dev.off()

```

### RDA of candidates

Probably don't include this as I really struggle with interpretation of its results.

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
Prepare R for RDA on laptop - Knowing Geography is a major factor we will try to discount it in the RDA
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
env.predictors <- env.predictors[, c('bio2', 'snow_days')]

# Load GenoType Data - choose 
gt_data <- read.PLINK("4-Downstream_Analyses/tables/fivepcs_topsnps.raw")
gt_data_tsv <- data.frame(as.matrix(gt_data))

# Coordinate data
samples_to_remove <- c("c_ll_ba_0216","c_ll_ba_0233","c_ll_cr_0211","h_ll_ba_0214","h_ll_ba_0215")
coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_wholeset.csv",
                          col_names = T, delim = ';')  %>% column_to_rownames(., var="id")
coord_table <- coord_table[!(row.names(coord_table) %in% samples_to_remove),]
coords <- data.frame(x=as.numeric(coord_table$longitude), y=as.numeric(coord_table$latitude),
                     row.names = rownames(coord_table))

env.predictors <- cbind(env.predictors, coords$x)
```
Run RDA
```{R}
rda <- rda(gt_data ~ ., data=env.predictors, scale=T)
rda <- rda(gt_data ~ . + Condition(coords$x) + Condition(coords$y), data=env.predictors, scale=T)
rda <- rda(gt_data ~ . + Condition(env.predictors$snow_days), data=env.predictors, scale=T)

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

Useful to add to a map
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

### Extract full list of genes

To extract the genes present in the candidate regions:
```{bash}
# On genomics-a
cd /home/ebazzicalupo/Selection_Eurasian_Lynx/Annotation

# bedtools to intersect bed with gff3 and get a bed of overlapping gene regions
bedtools intersect -wb \
 -a ../Intersect/total_intersect_candidate_fivepcs_windows.bed \
 -b <(awk -F"\t" '$3 == "gene" {print $0}' /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.gff3) |
 cut -f1,2,3,12 \
 > total_intersect_candidate_windows_fivepcs_intersect_genes.bed
# a total of 897 regions were found

# extract list of genes
cut -d'-' -f2 total_intersect_candidate_windows_intersect_genes.bed | 
 cut -d';' -f1 | sort -u \
 > total_intersect_candidate_windows_fivepcs_intersect_genes_list.txt
# a total of 782 gene IDs were found
# 104 of them have a LOC identifying number from the annotation
# -- CHECK what they mean -- #

# extract list of genes in the inversion
bedtools intersect -wb \
 -a ../Inversion/inversion.bed \
 -b <(awk -F"\t" '$3 == "gene" {print $0}' /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.gff3) |
 cut -f1,2,3,12 | cut -d'-' -f2 | cut -d';' -f1 | sort -u \
 > inversion_genes_list.txt
```
Download gene list to laptop
```{bash}
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/Annotation/total_intersect_candidate_windows_fivepcs_intersect_genes_list.txt Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/

scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/Annotation/inversion_genes_list.txt Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/
```
Input the list of gene IDs into Panther to make a comparison using the Felis catus annotation for GO-term representation.

Results are in 4-Downstream_Analyses/tables/PANTHER_Overrepresentation_Test_fivepcs_bioprocess_GO.txt

Genes in the inversion are: COP1 and PAPPA2

### Rank genes based on strongest signals

