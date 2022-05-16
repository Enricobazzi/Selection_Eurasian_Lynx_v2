---
title: "Prepare_Environmental_Data"
author: "Enrico"
date: "24/1/2022"
output: html_document
editor_options:
  chunk_output_type: console
---

In this markdown I will describe the process of extracting and selecting environmental data.

## Extracting

## Selecting environmental variables

### correlation and VIF between variables
Prepare R and load data
```{R}
library(tidyverse)
library("Hmisc")
library(corrplot)
library(viridis)
library(psych)
library(VIF)
library(fmsb)
library(usdm)

vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

varmatrix <- read_tsv("2-Prepare_Environmental_Data/WorldClim_table_persample.tsv") %>% column_to_rownames(., var="sample")
snowdata <- read_tsv("2-Prepare_Environmental_Data/Snow_table_persample.tsv") %>% column_to_rownames(., var="sample")
varmatrix2 <- as.matrix(cbind(varmatrix, snowdata))
varmatrix3 <- t(as.matrix(cbind(varmatrix, snowdata)))
```
Extract Matrix with correlations
```{R}
pdf(file = paste0("2-Prepare_Environmental_Data/plots/allvars_pairs.panels.pdf"),
    width = 13,
    height = 8)
pairs.panels(varmatrix2, scale=T, cex.cor = 1.5)
dev.off()
```
Selection of 4 main categories based on correlation values (<0.7):
 - 1,2,3,4,6,7,9,11,12,14,15,17,19
 - 5,10
 - 8
 - 13,16,18
 - Snow days
 - Jan depth
Selecting 5, 8 and 13 and snow variables for easy interpretation, discarding 3 and 4 for harder interpretation, and 12 for correlation with 13, we need to choose between:
 - 1,2,4,6,7,9,11,14,15,17,19

Reitarative VIF method
```{R}
# all variables
varmatrix_short <- varmatrix2[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
# all competing variables (removed 3,4,10,12,13,16,18 for reasons mentioned above)
varmatrix_short <- varmatrix2[,c(1,2,7,6,5,8,9,11,13,14,15,17,20,21,3,4)]

# manually do by removing highest VIF from c() or use function below (depends on order of columns) 
varmatrix_short <- varmatrix2[,c(1,2,7,6,5,8,9,11,13,14,15,17,20,21)]
usdm::vifcor(varmatrix_short, th=10)

# final set of variables - loose
vif_func(varmatrix_short, trace =T, thresh = 4)
varmatrix_final <- varmatrix2[,c(2,5,6,8,13,20,21)]

usdm::vifcor(varmatrix_final, th=10)
#   Variables   VIF
#      bio2     3.976499
#      bio5     2.158984
#      bio6     3.750249
#      bio8     1.594988
#     bio13     1.412701
# jan_depth     2.236679
# snow_days     3.329262


# final set of variables - strict
vif_func(varmatrix_short, trace =T, thresh = 2)
varmatrix_final <- varmatrix2[,c(5,6,8,13,20)]

usdm::vifcor(varmatrix_final, th=10)
# Variables      VIF
#      bio6 1.261347
#      bio5 1.467656
#      bio8 1.328688
#     bio13 1.218006
# jan_depth 1.241692

# Export table of selected variables only
export_final <- data.frame(varmatrix_final) %>% rownames_to_column(., var = "sample")

write_tsv(export_final, "2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv") 
```
Correlation plot of uncorrelated set of variables
```{R}
varmatrix <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv") %>% column_to_rownames(., var="sample")
pdf(file = paste0("2-Prepare_Environmental_Data/plots/uncorrelated_variables_pairs.panels.pdf"),
    width = 13,
    height = 8)
pairs.panels(varmatrix, scale=T, cex.cor = 1.5)
dev.off()
```
Correlation of variables across the whole distributional range
```{R}
library(raster)
library(rgdal)
library(grDevices)
library(sf)
library(viridis)
library(RColorBrewer)
library(dismo)
library(tidyverse)
library("Hmisc")
library(corrplot)
library(viridis)
library(psych)
library(VIF)
library(fmsb)
library(usdm)


# Create raster with all variables
worldclim <- getData("worldclim", var = "bio", res = 10)
jan_mean_depth <- raster("2-Prepare_Environmental_Data/tables/jan_mean_depth.tif")
jan_mean_depth_new <- projectRaster(jan_mean_depth, worldclim)
mean_snow_days <- raster("2-Prepare_Environmental_Data/tables/mean_snow_days.tif")
mean_snow_new <- projectRaster(mean_snow_days, worldclim)

clim.layer <- stack(worldclim,jan_mean_depth_new,mean_snow_new)

# crop for distributional range only - depe
distr.map <- readOGR("~/Downloads/redlist_species_data_1f4a1a8f-31fa-48ee-8678-7435f90a8ff9/data_0.shp")
distr.map <- readOGR("~/Downloads/redlist_species_data_5a932080-5312-4e23-b0d3-929acfb17324/data_0.shp")

r2 <- crop(clim.layer, extent(distr.map))
r3 <- mask(r2, distr.map)
ext <- extent(r3)

#To create 1000 random points
set.seed(123)
backgr <- randomPoints(r3, 1000, ext=ext)
backgrvals <- data.frame(extract(r3, backgr))

#To plot it you can use
plot(r3[[1]])
points(backgr, col='black')

# Correlation plots
pdf(file = paste0("2-Prepare_Environmental_Data/plots/all_variables_pairs.panels.pdf"),
    width = 13,
    height = 8)
pairs.panels(backgrvals, scale=T, cex.cor = 1.5)
dev.off()
# The following "groups" of r>0.7
# 1,4,6,7,9,11
# 2,14,15,17,19(not with 2)
# 5,10
# 8
# 13,16,18
# 3 -> only with 4 - weird being with only one from one group but not rest
# 12 -> with 13,14,16,17 + 0.69 with 18,19 - weird being with more than 1 group
# snow variables highly correlated across range
# of these groups our sample data VIF selection has selected one candidate from each
# except for bio3

usdm::vifcor(backgrvals, th=10)
varmatrix_short <- as.matrix(backgrvals[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)])

# final set of variables - loose
vif_func(varmatrix_short, trace=T, thresh = 4)
```

## PCA of samples based on variables

```{R}
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(viridis)
library(RColorBrewer)

varmatrix <- read_tsv("2-Prepare_Environmental_Data/WorldClim_table_persample.tsv") %>%
  column_to_rownames(., var="sample")
snowdata <- read_tsv("2-Prepare_Environmental_Data/Snow_table_persample.tsv") %>%
  column_to_rownames(., var="sample")
var.data <- cbind(varmatrix, snowdata)

row.names(var.data)
samples_to_remove <- c("c_ll_ba_0216","c_ll_ba_0233","c_ll_cr_0211","h_ll_ba_0214","h_ll_ba_0215","c_ll_no_0065", "c_ll_no_0075",
                       "c_ll_no_0076", "c_ll_no_0077", "c_ll_no_0078", "c_ll_no_0079", "c_ll_no_0080", "c_ll_no_0081", "c_ll_no_0082",
                       "c_ll_vl_0137", "c_ll_og_0181", "c_ll_og_0187", "c_ll_cr_0205", "c_ll_cr_0206", 
                       "c_ll_cr_0207", "c_ll_cr_0208", "c_ll_cr_0209", "c_ll_cr_0212", "c_ll_ba_0224", "c_ll_ba_0226",
                       "c_ll_ba_0227", "c_ll_ba_0228", "c_ll_ba_0229", "c_ll_ba_0230", "c_ll_tu_0154",
                       "c_ll_po_0001", "c_ll_po_0002", "c_ll_po_0003", "c_ll_po_0011", "c_ll_po_0014", "c_ll_po_0019",
                       "c_ll_po_0105", "c_ll_po_0106", "c_ll_po_0150")

var.data <- var.data[!(row.names(var.data) %in% samples_to_remove),]

loc <- rep(NA, NROW(var.data))
loc[grep("ba", rownames(var.data))] <- "Balkans"
loc[grep("ca", rownames(var.data))] <- "Caucasus"
loc[grep("cr", rownames(var.data))] <- "Carpathians"
loc[grep("ka", rownames(var.data))] <- "Mongolia"
loc[grep("ki", rownames(var.data))] <- "Kirov"
loc[grep("la", rownames(var.data))] <- "Latvia"
loc[grep("no", rownames(var.data))] <- "Norway"
loc[grep("po", rownames(var.data))] <- "NE-Poland"
loc[grep("og", rownames(var.data))] <- "Mongolia"
loc[grep("to", rownames(var.data))] <- "Mongolia"
loc[grep("tu", rownames(var.data))] <- "Tuva"
loc[grep("ur", rownames(var.data))] <- "Urals"
loc[grep("vl", rownames(var.data))] <- "Vladivostok"
loc[grep("ya", rownames(var.data))] <- "Yakutia"

num <- rep(NA, NROW(var.data))
num[grep("ba", rownames(var.data))] <- 1
num[grep("ca", rownames(var.data))] <- 3
num[grep("cr", rownames(var.data))] <- 2
num[grep("ka", rownames(var.data))] <- 6
num[grep("ki", rownames(var.data))] <- 4
num[grep("la", rownames(var.data))] <- 5
num[grep("no", rownames(var.data))] <- 8
num[grep("po", rownames(var.data))] <- 7
num[grep("og", rownames(var.data))] <- 6
num[grep("to", rownames(var.data))] <- 6
num[grep("tu", rownames(var.data))] <- 9
num[grep("ur", rownames(var.data))] <- 10
num[grep("vl", rownames(var.data))] <- 11
num[grep("ya", rownames(var.data))] <- 12

ola <- data.frame(sample = rownames(var.data), pop = loc, n = num)

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


pca <- PCA(var.data, scale.unit = TRUE, ncp = 5, graph = TRUE)

# Screeplot
pdf(file = paste0("2-Prepare_Environmental_Data/plots/screeplot_allvars_pca.pdf"),
    width = 8,
    height = 8)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 60))
dev.off()

# cos2 (quality of representation) of each variable
# if only 2 PCs (~70% of variance)
fviz_cos2(pca, choice = "var", axes = 1:2)
# if only 3 PCs (~80% of variance)
fviz_cos2(pca, choice = "var", axes = 1:4)
# if 5 PCs (~95% of variance)
fviz_cos2(pca, choice = "var", axes = 1:5)

for(n in 2:5){
  pdf(file = paste0("2-Prepare_Environmental_Data/plots/cos2_pca_pcs1to", n,".pdf"),
    width = 8,
    height = 8)
  fviz_cos2(pca, choice = "var", axes = 1:n)
  dev.off()
}

# Color by cos2 values: quality on the factor map
fviz_pca_var(pca, col.var = "cos2", axes = 2:3, 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
             )
var <- get_pca_var(pca)
corrplot(var$contrib, is.corr=FALSE)

# Contributions of variables to PC1
fviz_contrib(pca, choice = "var", axes = 1, top = 21)
# Contributions of variables to PC2
fviz_contrib(pca, choice = "var", axes = 2, top = 21)
# Contributions of variables to PC3
fviz_contrib(pca, choice = "var", axes = 3, top = 21)
# Contributions of variables to PC4
fviz_contrib(pca, choice = "var", axes = 4, top = 21)
# Contributions of variables to PC5
fviz_contrib(pca, choice = "var", axes = 5, top = 21)

desc <- dimdesc(pca, axes = c(1,2,3,4,5), proba = 0.05)
desc$Dim.1
desc$Dim.2
desc$Dim.3
desc$Dim.4
desc$Dim.5

ind <- get_pca_ind(pca)

# EXTRACT TABLE OF PC VALUES of each sample
twopcs.dataframe <- data.frame(PC1=ind$coord[,1], PC2=ind$coord[,2]) %>%
                     rownames_to_column("sample")

write_tsv(twopcs.dataframe, "2-Prepare_Environmental_Data/twopcs_data_matrix.tsv")


fivepcs.dataframe <- data.frame(PC1=ind$coord[,1], PC2=ind$coord[,2], PC3=ind$coord[,3], 
                                PC4=ind$coord[,4], PC5=ind$coord[,5]) %>%
                     rownames_to_column("sample")

write_tsv(fivepcs.dataframe, "2-Prepare_Environmental_Data/fivepcs_data_matrix.tsv")

# EXTRACT TABLE OF PC VALUES for each population
# change the order of the columns so that they are in alphabetical order 
# (as in allele counts file: ca - ki - la - mo - tu - ur - vl - ya)
populations <- c("ca", "ki", "la", "mo", "tu", "ur", "vl", "ya")
final.table <- data.frame(variable=c("PC1", "PC2", "PC3", "PC4", "PC5"))
for(p in 1:length(populations)){
  # population
  POP <- populations[p]
  if (POP=="mo"){
    pop.sample.n <- grep("ka|og|to", rownames(ind$coord))
  } else {
  # rows of population
  pop.sample.n <- grep(POP, rownames(ind$coord))
  }
  # create a table just for the population adding samples one by one
  pop.table <- data.frame()
  for(n in 1:length(pop.sample.n)){
    N <- pop.sample.n[n]
    row <- data.frame(PC1=ind$coord[N,1], PC2=ind$coord[N,2], PC3=ind$coord[N,3], 
                                PC4=ind$coord[N,4], PC5=ind$coord[N,5])
    pop.table <- rbind(pop.table, row)
  }
  # get mean value of population for each PC population's v
  pc1.mean <- mean(pop.table$PC1)
  pc2.mean <- mean(pop.table$PC2)
  pc3.mean <- mean(pop.table$PC3)
  pc4.mean <- mean(pop.table$PC4)
  pc5.mean <- mean(pop.table$PC5)
  pop.column.vals <- c(pc1.mean, pc2.mean, pc3.mean, pc4.mean, pc5.mean)
  final.table <- cbind(final.table, data.frame(V=pop.column.vals))
}
colnames(final.table) <- c("variable", populations)

write_tsv(final.table, "2-Prepare_Environmental_Data/fivepcs_populations_data_matrix.tsv")
```
