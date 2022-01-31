library(raster)
library(rgdal)
library(grDevices)
library(sf)
library(viridis)
library(RColorBrewer)

samples_to_remove <- c("c_ll_ba_0216","c_ll_ba_0233","c_ll_cr_0211","h_ll_ba_0214","h_ll_ba_0215")

clim.data <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv", col_names = T) %>%
  column_to_rownames(., var="sample")

coord_table <- read_delim("~/Dropbox/LL_LC_LR_Databases/LL_coords/csv_LL_selection_coords_wholeset.csv",
                          col_names = T, delim = ';')  %>% column_to_rownames(., var="id")
coord_table <- coord_table[!(row.names(coord_table) %in% samples_to_remove),]

clim.points <- data.frame(clim.data, x=as.numeric(coord_table$longitude), y=as.numeric(coord_table$latitude)) %>% 
  rownames_to_column(., var="ID")

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

worldclim <- getData("worldclim", var = "bio", res = 10)

data <- worldclim[[2]]

binary_data <- data > 0

extent <- c(-10, 180, 20, 90) 

data.crop <- crop(binary_data, extent)
distr.map <- readOGR("~/Downloads/redlist_species_data_1f4a1a8f-31fa-48ee-8678-7435f90a8ff9/data_0.shp")
r2 <- crop(data.crop, extent(distr.map))
r3 <- mask(r2, distr.map)

cuts=c(0.99,1,1.01) #set breaks
pal <- colorRampPalette(c("black","black"))
pal2 <- colorRampPalette(c("darkgreen","darkgreen"))

pdf("2-Prepare_Environmental_Data/plots/distribution_map.pdf", width = 14)
plot(data.crop, alpha=0.4, breaks=cuts, col = pal(3), integrate=T)
plot(r3, breaks=cuts, col = pal2(3), integrate=T, add=T)
plot(distr.map, lwd=1.5, add=T)
points(clim.points$x,clim.points$y, pch=21, lwd=1.5, cex=1.2, col="black", bg=bg[num]) # the lynxes
dev.off()
