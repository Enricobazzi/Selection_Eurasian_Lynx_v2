# An inversion BED was created manually like this:
# echo -E "scaffold_17_arrow_ctg1\t17635000\t18355000" > Inversion/inversion.bed

# # Then I ran PCA with PLINK using all samples like this:
# cd /home/ebazzicalupo/Selection_Eurasian_Lynx/
# bedtools intersect -header -a VCF/ll_wholegenome_LyCa_ref.sorted.filter7.vcf -b Inversion/inversion.bed \
#  > Inversion/inversion.vcf

# # I manually created a file listing samples to remove for plink with following format:
# # c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# # c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# # c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# # h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# # h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

# plink_1.9 --vcf VCF/inversion.vcf \
#  --double-id --allow-extra-chr --set-missing-var-ids @:# \
#  --remove VCF/samplestoremove.txt \
#  --pca --out Inversion/inversion_pca

# Donwnload to laptop
# scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/Selection_Eurasian_Lynx/Inversion/inversion_pca.eigen\* \
# Documents/Selection_Eurasian_Lynx_v2/4-Downstream_Analyses/tables/

library(tidyverse)
library(RColorBrewer)
library(viridis)

# import eigen vec and val
pca <- read_table2("4-Downstream_Analyses/tables/superoutlier_persample.eigenvec",
                   col_names = FALSE)
eigenval <- scan("4-Downstream_Analyses/tables/superoutlier_persample.eigenval")

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
pdf(file = paste0("4-Downstream_Analyses/plots/inversion_pca1_pca2.pdf"),
    width = 8,
    height = 8)

ggplot() + 
  geom_point(data=pca, aes(PC1, PC2), size = 4, pch=21, bg=loc, col="gray32") +
  coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  scale_color_manual(values=cols) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

dev.off()

### PIE CHARTS of each population's allele frequency of inversion

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
