neutral_pca.rast.crop
candidate_pca.rast.crop

extent <- c(-10, 180, 20, 80)

neutral_pca.rast.crop <- crop(neutral_pca.rast, extent)
candidate_pca.rast.crop <- crop(candidate_pca.rast, extent)

neutral_clim.rast <- na.omit(getValues(neutral_clim.trans))
candidate_clim.rast <- na.omit(getValues(candidate_clim.trans))

env_trns <- data.frame(getValues(neutral_clim.trans)) %>% rownames_to_column("cells") %>% drop_na()

####

neutral_pca <- prcomp(neutral_clim.rast)
candidate_pca <- prcomp(candidate_clim.rast)

candidate_pca.rast <- predict(candidate_clim.trans, candidate_pca, index=1:5)
neutral_pca.rast <- predict(neutral_clim.trans, neutral_pca, index=1:5)
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

####
# NEUTRAL PCA LEGEND
neutral_legend.colors <- na.omit(getValues(neutral_pca.rast[[1:3]]))
neutral_colors <- c()
for (r in 1:NROW(neutral_legend.colors)){
  print(r)
  red <- neutral_legend.colors[r,1]/255
  green <- neutral_legend.colors[r,2]/255
  blue <- neutral_legend.colors[r,3]/255
  colo <- rgb(red,green,blue,alpha=1)
  neutral_colors <- c(neutral_colors,colo)
}

sample(1:100, 3, replace=TRUE)

plotcoords=data.frame(neutral_pca$x)
envloading1=data.frame(neutral_pca$rotation)
envloading <- envloading1

# prettier plot
PCA <- ggplot()+
  #scale_color_viridis(discrete = soils$Core)+
  #scale_color_manual(values = wes_palette("BottleRocket2", 100, type="continuous"))+
  geom_point(aes(x=plotcoords$PC1, y=plotcoords$PC2), size = 4, shape=20, color=neutral_colors)+
  #scale_color_manual(values = colors)+
  scale_x_continuous(limits=c(-0.25,1.15))+
  geom_segment(aes(x=c(0,0,0,0,0,0,0,0,0,0,0), y=c(0,0,0,0,0,0,0,0,0,0,0), xend= c(envloading[1,1]/2,envloading[2,1]/2,envloading[3,1]/2,envloading[4,1]/2,envloading[5,1]/2,envloading[6,1]/2,envloading[7,1]/2,envloading[8,1]/2,envloading[9,1]/2,envloading[10,1]/2,envloading[11,1]/2), yend=c(envloading[1,2]/2,envloading[2,2]/2,envloading[3,2]/2,envloading[4,2]/2,envloading[5,2]/2,envloading[6,2]/2,envloading[7,2]/2,envloading[8,2]/2,envloading[9,2]/2,envloading[10,2]/2,envloading[11,2]/2)), color="black", arrow=arrow(angle = 20, length = unit(0.4,"cm"), ends = "last", type = "open"), size = 1.2)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("PC 1 -- 63% of variance")+
  ylab("PC 2 -- 18% of variance")+
  #labs(colour = "Core in Transect")+
  # geom_text(aes(x=plotcoords$PC1, y=plotcoords$PC2, label = pcadata$Sample, color = pcadata$Transect), nudge_x =0.25)+
  theme_void() + # remove background, grid, numeric labels
  annotate("text", 
           x=c((envloading[1,1]+0.3)/2,
               (envloading[2,1]+0.55)/2,
               (envloading[5,1]-0.1)/2,
               (envloading[10,1]+0.35)/2),
           y=c((envloading[1,2]-0.08)/2,
               envloading[2,2]/2,
               (envloading[5,2]-0.08)/2,
               (envloading[10,2]-0.03)/2),
           label =  c("x-coord","y-coord","bio3","bio15"), color = "black", size = 9)

PCA

pdf(file = paste0("4-Downstream_Analyses/plots/GDM_neutral_legend.pdf"),
    width = 4,
    height = 4)
PCA
dev.off()

# CANDIDATES PCA LEGEND
candidate_legend.colors <- na.omit(getValues(candidate_pca.rast[[1:3]]))
candidate_colors <- c()
for (r in 1:NROW(candidate_legend.colors)){
   print(r)
   red <- candidate_legend.colors[r,1]/255
   green <- candidate_legend.colors[r,2]/255
   blue <- candidate_legend.colors[r,3]/255
   colo <- rgb(red,green,blue,alpha=1)
   candidate_colors <- c(candidate_colors,colo)
}

sample(1:100, 3, replace=TRUE)

plotcoords=data.frame(candidate_pca$x)
envloading1=data.frame(candidate_pca$rotation)
envloading <- envloading1

# prettier plot
PCA <- ggplot()+
  #scale_color_viridis(discrete = soils$Core)+
  #scale_color_manual(values = wes_palette("BottleRocket2", 100, type="continuous"))+
  geom_point(aes(x=plotcoords$PC1, y=plotcoords$PC2), size = 4, shape=20, color=candidate_colors)+
  #scale_color_manual(values = colors)+
  scale_x_continuous(limits=c(-0.25,0.7))+
  geom_segment(aes(x=c(0,0,0,0,0,0,0,0,0), y=c(0,0,0,0,0,0,0,0,0), xend= c(envloading[1,1]/2,envloading[2,1]/2,envloading[3,1]/2,envloading[4,1]/2,envloading[5,1]/2,envloading[6,1]/2,envloading[7,1]/2,envloading[8,1]/2,envloading[9,1]/2), yend=c(envloading[1,2]/2,envloading[2,2]/2,envloading[3,2]/2,envloading[4,2]/2,envloading[5,2]/2,envloading[6,2]/2,envloading[7,2]/2,envloading[8,2]/2,envloading[9,2]/2)), color="black", arrow=arrow(angle = 20, length = unit(0.4,"cm"), ends = "last", type = "open"), size = 1.2)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("PC 1 -- 38% of variance")+
  ylab("PC 2 -- 18% of variance")+
  #labs(colour = "Core in Transect")+
  # geom_text(aes(x=plotcoords$PC1, y=plotcoords$PC2, label = pcadata$Sample, color = pcadata$Transect), nudge_x =0.25)+
  theme_void() + # remove background, grid, numeric labels
  annotate("text", 
           x=c((envloading[1,1]+0.1)/2,
               (envloading[2,1]+0.42)/2,
               (envloading[4,1]+0.18)/2,
               (envloading[5,1]-0.14)/2,
               (envloading[9,1]+0.72)/2), 
           y=c((envloading[1,2]-0.08)/2,
               envloading[2,2]/2,
               (envloading[4,2]-0.09)/2,
               (envloading[5,2]-0.09)/2,
               (envloading[9,2]-0.1)/2), 
           label =  c("x-coord","y-coord","bio2","bio3", "mean snow days"), color = "black", size = 9)

#PCA

pdf(file = paste0("4-Downstream_Analyses/plots/GDM_candidate_legend.pdf"),
    width = 4,
    height = 4)
PCA
dev.off()

######



PCA <- ggplot()+
  #scale_color_viridis(discrete = soils$Core)+
  #scale_color_manual(values = wes_palette("BottleRocket2", 100, type="continuous"))+
  geom_point(aes(x=plotcoords$PC1, y=plotcoords$PC3), size = 4, shape=20, color=candidate_colors)+
  #scale_color_manual(values = colors)+
  scale_x_continuous(limits=c(-0.25,0.7))+
  geom_segment(aes(x=c(0,0,0,0,0,0,0,0,0), y=c(0,0,0,0,0,0,0,0,0), xend= c(envloading[1,1]/2,envloading[2,1]/2,envloading[3,1]/2,envloading[4,1]/2,envloading[5,1]/2,envloading[6,1]/2,envloading[7,1]/2,envloading[8,1]/2,envloading[9,1]/2), yend=c(envloading[1,3]/2,envloading[2,3]/2,envloading[3,3]/2,envloading[4,3]/2,envloading[5,3]/2,envloading[6,3]/2,envloading[7,3]/2,envloading[8,3]/2,envloading[9,3]/2)), color="black", arrow=arrow(angle = 20, length = unit(0.4,"cm"), ends = "last", type = "open"), size = 1.2)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("PC 1 -- 38% of variance")+
  ylab("PC 3 -- ??% of variance")+
  #labs(colour = "Core in Transect")+
  # geom_text(aes(x=plotcoords$PC1, y=plotcoords$PC2, label = pcadata$Sample, color = pcadata$Transect), nudge_x =0.25)+
  theme_void() + # remove background, grid, numeric labels
  annotate("text", 
           x=c((envloading[1,1])/2,
               (envloading[2,1])/2,
               (envloading[4,1])/2,
               (envloading[5,1])/2,
               (envloading[9,1])/2), 
           y=c((envloading[1,3])/2,
               envloading[2,3]/2,
               (envloading[4,3])/2,
               (envloading[5,3])/2,
               (envloading[9,3])/2), 
           label =  c("x-coord","y-coord","bio2","bio3", "mean snow days"), color = "black", size = 9)

PCA
