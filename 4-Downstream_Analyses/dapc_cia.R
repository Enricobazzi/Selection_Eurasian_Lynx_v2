# DAPC of neutral
# kmeans
grp <- find.clusters(adegen_neutral_raw, n.pca = 103, n.clust=4, max.n.clust=20)
neutral_grps <- grp$grp
# cross validate NPCs to use in DAPC
xval <- xvalDapc(adegen_neutral_raw, neutral_grps, n.pca.max = 200,
                 result = "overall",
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)
# DAPC
neutral_DAPC <- dapc(adegen_neutral_raw, neutral_grps, n.pca = xval$DAPC$n.pca, n.da = xval$DAPC$n.da)
# plot DAPC
myCol <- c(brewer.pal(n = 7, name = "Dark2"))
scatter(neutral_DAPC, col=myCol, scree.da=FALSE, xax = 1, yax = 2,
        cell=1.5, cex=2, bg="white",cstar=1)
scatter(neutral_DAPC, col=myCol, scree.da=FALSE, xax = 1, yax = 3,
        cell=1.5, cex=2, bg="white",cstar=1)
# get matrix of DAPC
neutral_dapc_matrix <- neutral_DAPC$ind.coord


# DAPC of candidates
grp <- find.clusters(adegen_candidate_raw, n.pca = 103, n.clust=6, max.n.clust = 20)

candidate_grps <- grp$grp

xval <- xvalDapc(adegen_candidate_raw, candidate_grps, n.pca.max = 200,
                 result = "overall",
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)

candidate_DAPC <- dapc(adegen_candidate_raw, candidate_grps, n.pca = xval$DAPC$n.pca, n.da = xval$DAPC$n.da)
candidate_DAPC$posterior

myCol <- c(brewer.pal(n = 8, name = "Dark2"))
scatter(candidate_DAPC, col=myCol, scree.da=FALSE, xax = 1, yax = 2,
        cell=1.5, cex=2, bg="white",cstar=1)

candidate_dapc_matrix <- candidate_DAPC$ind.coord

# get snps responsible for axis:
set.seed(4)
contrib <- loadingplot(candidate_DAPC$var.contr, axis=1,
                       thres=.008, lab.jitter=1)

# made4
cia1 <- cia(t(neutral_dapc_matrix), t(candidate_dapc_matrix), cia.nf=3, cia.scan=FALSE, nsc=TRUE)
cia1$coinertia$RV
s.match(cia1$coinertia$mX, cia1$coinertia$mY, clabel = 0.2)
summary(cia1$coinertia)
plot.coinertia
dapc_cia <- data.frame(sample=(rownames(cia1$coinertia$mX)), 
                       start_x=(cia1$coinertia$mX$NorS1), start_y=(cia1$coinertia$mX$NorS2),
                       end_x=(cia1$coinertia$mY$NorS1), end_y=(cia1$coinertia$mY$NorS2),
                       neutral_grp=(as.numeric(neutral_grps)+20), candidate_grp=(as.numeric(candidate_grps)+20))

loc <- rep(NA, length(dapc_cia$sample))
loc[grep("ca", dapc_cia$sample)] <- "#B8860b"
loc[grep("po", dapc_cia$sample)] <- viridis_pal()(5)[3]
loc[grep("ba", dapc_cia$sample)] <- "#A035AF"
loc[grep("cr", dapc_cia$sample)] <- brewer.pal(12,"Paired")[9]
loc[grep("no", dapc_cia$sample)] <- viridis_pal()(5)[2]
loc[grep("ka", dapc_cia$sample)] <- brewer.pal(12,"Paired")[7]
loc[grep("ki", dapc_cia$sample)] <- viridis_pal()(5)[1]
loc[grep("la", dapc_cia$sample)] <- brewer.pal(12,"Paired")[3]
loc[grep("to", dapc_cia$sample)] <- brewer.pal(12,"Paired")[7]
loc[grep("og", dapc_cia$sample)] <- brewer.pal(12,"Paired")[7]
loc[grep("tu", dapc_cia$sample)] <- brewer.pal(12,"Paired")[8]
loc[grep("ur", dapc_cia$sample)] <- "#0F4909"
loc[grep("vl", dapc_cia$sample)] <- brewer.pal(12,"Paired")[5]
loc[grep("ya", dapc_cia$sample)] <- brewer.pal(12,"Paired")[6]

dapc_cia <- data.frame(dapc_cia, color=(loc))

ggplot() + 
  geom_point(data=dapc_cia, mapping=aes(x=start_x, y=start_y), size=3, fill=loc) +
  geom_segment(data=dapc_cia, mapping=aes(x=start_x, y=start_y, xend=end_x, yend=end_y), size=0.4, color="darkgrey") + 
  geom_point(data=dapc_cia, mapping=aes(x=end_x, y=end_y), size=5, shape=21, fill=loc) +
  theme_minimal()

igraph::compare(candidate_grps,neutral_grps, method="vi")
