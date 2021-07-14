###################################
#### Hierarchy metrics ####
###################################


library(EloRating)
library(steepness)
library(compete)
library(ggplot2)
library(GGally)

load('./Data/hierarchyMetricsData.RData')

dominance<-lapply(data,convert.dom)

sge.I.hierarchy.metrics <- matrix(nrow = 18, ncol = 6)
count <- 1
for (i in unique(names(dominance))){
  tmp <- as.data.frame(dominance[i])
  tmp <- tmp[order(tmp[,1]),]
  tmp <- subset(tmp, tmp[,3] != tmp[,4])
  seqcheck(tmp[,3],tmp[,4],tmp[,1])
  elo<-elo.seq(tmp[,3],tmp[,4],tmp[,1])
  elo.mat <- creatematrix(elo)
  a <- extract_elo(elo)
  b <- steeptest(elo.mat, rep = 1000)$Stp
  c <- devries(elo.mat)$`h-modified`
  d <- dci(elo.mat)
  e <- phi(elo.mat)
  f <- ttri(elo.mat)$ttri
  g <- stab_elo(elo)
  sge.I.hierarchy.metrics[count,] <- c(b,c,d,e,f,g)
  count <- count + 1
}
sge.I.hierarchy.metrics <- as.data.frame(sge.I.hierarchy.metrics)
sge.I.hierarchy.metrics$phase <- rep(c(1,2), times = c(9,9))
rownames(sge.I.hierarchy.metrics) <- names(dominance)
colnames(sge.I.hierarchy.metrics) <- c('steepness','linearity','DC','phi','tt','stability','phase')
sge.I.hierarchy.metrics$group <- rownames(sge.I.hierarchy.metrics)

# plot correlations between elos, DS, and ordinal ranks for SGE I
ggpairs(as.data.frame(sge.I.hierarchy.metrics[,c(1,3:4,6)]), title = 'Phases 1 + 2')
ggpairs(as.data.frame(sge.I.hierarchy.metrics[1:9,c(1,3:4,6)]), title = 'Phase 1')
ggpairs(as.data.frame(sge.I.hierarchy.metrics[10:18,c(1,3:4,6)]), title = 'Phase 2')
