####################
#### fsgea code ####
####################

## load libraries
library(msigdbr)
library(fgsea)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(tidyverse)

## load data and pathways
pathways.hallmark <- gmtPathways("./Data/h.all.v7.2.symbols.gmt")
load('./Data/fgseaData.RData')

## preprocess data for fgsea function
CD4_agAsym_GSEAgenelist <- gc.agg.Diff_CD4
CD8_agAsym_GSEAgenelist <- gc.agg.Diff_CD8
CD20_agAsym_GSEAgenelist <- gc.agg.Diff_CD20
CD16_agAsym_GSEAgenelist <- gc.agg.Diff_CD16
CD14_agAsym_GSEAgenelist <- gc.agg.Diff_CD14
CD4_elo_GSEAgenelist <- elo_CD4
CD8_elo_GSEAgenelist <- elo_CD8
CD20_elo_GSEAgenelist <- elo_CD20
CD16_elo_GSEAgenelist <- elo_CD16
CD14_elo_GSEAgenelist <- elo_CD14

CD4_agAsym_GSEAgenelist$stBeta     <- CD4_agAsym_GSEAgenelist$beta_behaveVar / sqrt(CD4_agAsym_GSEAgenelist$var_behaveVar)
CD4_agAsym_GSEAgenelist            <- CD4_agAsym_GSEAgenelist[order(CD4_agAsym_GSEAgenelist$stBeta),]
CD4_agAsym_GSEAgenelist$geneSymbol <- rownames(CD4_agAsym_GSEAgenelist)
CD4_agAsym_GSEAgenelist2           <- as.vector(CD4_agAsym_GSEAgenelist$stBeta)
names(CD4_agAsym_GSEAgenelist2)    <- CD4_agAsym_GSEAgenelist$geneSymbol

CD4_elo_GSEAgenelist$stBeta     <- CD4_elo_GSEAgenelist$beta_behaveVar / sqrt(CD4_elo_GSEAgenelist$var_behaveVar)
CD4_elo_GSEAgenelist            <- CD4_elo_GSEAgenelist[order(CD4_elo_GSEAgenelist$stBeta),]
CD4_elo_GSEAgenelist$geneSymbol <- rownames(CD4_elo_GSEAgenelist)
CD4_elo_GSEAgenelist2           <- as.vector(CD4_elo_GSEAgenelist$stBeta)
names(CD4_elo_GSEAgenelist2)    <- CD4_elo_GSEAgenelist$geneSymbol

CD8_agAsym_GSEAgenelist$stBeta     <- CD8_agAsym_GSEAgenelist$beta_behaveVar / sqrt(CD8_agAsym_GSEAgenelist$var_behaveVar)
CD8_agAsym_GSEAgenelist            <- CD8_agAsym_GSEAgenelist[order(CD8_agAsym_GSEAgenelist$stBeta),]
CD8_agAsym_GSEAgenelist$geneSymbol <- rownames(CD8_agAsym_GSEAgenelist)
CD8_agAsym_GSEAgenelist2           <- as.vector(CD8_agAsym_GSEAgenelist$stBeta)
names(CD8_agAsym_GSEAgenelist2)    <- CD8_agAsym_GSEAgenelist$geneSymbol

CD8_elo_GSEAgenelist$stBeta     <- CD8_elo_GSEAgenelist$beta_behaveVar / sqrt(CD8_elo_GSEAgenelist$var_behaveVar)
CD8_elo_GSEAgenelist            <- CD8_elo_GSEAgenelist[order(CD8_elo_GSEAgenelist$stBeta),]
CD8_elo_GSEAgenelist$geneSymbol <- rownames(CD8_elo_GSEAgenelist)
CD8_elo_GSEAgenelist2           <- as.vector(CD8_elo_GSEAgenelist$stBeta)
names(CD8_elo_GSEAgenelist2)    <- CD8_elo_GSEAgenelist$geneSymbol

CD16_agAsym_GSEAgenelist$stBeta     <- CD16_agAsym_GSEAgenelist$beta_behaveVar / sqrt(CD16_agAsym_GSEAgenelist$var_behaveVar)
CD16_agAsym_GSEAgenelist            <- CD16_agAsym_GSEAgenelist[order(CD16_agAsym_GSEAgenelist$stBeta),]
CD16_agAsym_GSEAgenelist$geneSymbol <- rownames(CD16_agAsym_GSEAgenelist)
CD16_agAsym_GSEAgenelist2           <- as.vector(CD16_agAsym_GSEAgenelist$stBeta)
names(CD16_agAsym_GSEAgenelist2)    <- CD16_agAsym_GSEAgenelist$geneSymbol

CD16_elo_GSEAgenelist$stBeta     <- CD16_elo_GSEAgenelist$beta_behaveVar / sqrt(CD16_elo_GSEAgenelist$var_behaveVar)
CD16_elo_GSEAgenelist            <- CD16_elo_GSEAgenelist[order(CD16_elo_GSEAgenelist$stBeta),]
CD16_elo_GSEAgenelist$geneSymbol <- rownames(CD16_elo_GSEAgenelist)
CD16_elo_GSEAgenelist2           <- as.vector(CD16_elo_GSEAgenelist$stBeta)
names(CD16_elo_GSEAgenelist2)    <- CD16_elo_GSEAgenelist$geneSymbol

CD20_agAsym_GSEAgenelist$stBeta     <- CD20_agAsym_GSEAgenelist$beta_behaveVar / sqrt(CD20_agAsym_GSEAgenelist$var_behaveVar)
CD20_agAsym_GSEAgenelist            <- CD20_agAsym_GSEAgenelist[order(CD20_agAsym_GSEAgenelist$stBeta),]
CD20_agAsym_GSEAgenelist$geneSymbol <- rownames(CD20_agAsym_GSEAgenelist)
CD20_agAsym_GSEAgenelist2           <- as.vector(CD20_agAsym_GSEAgenelist$stBeta)
names(CD20_agAsym_GSEAgenelist2)    <- CD20_agAsym_GSEAgenelist$geneSymbol

CD20_elo_GSEAgenelist$stBeta     <- CD20_elo_GSEAgenelist$beta_behaveVar / sqrt(CD20_elo_GSEAgenelist$var_behaveVar)
CD20_elo_GSEAgenelist            <- CD20_elo_GSEAgenelist[order(CD20_elo_GSEAgenelist$stBeta),]
CD20_elo_GSEAgenelist$geneSymbol <- rownames(CD20_elo_GSEAgenelist)
CD20_elo_GSEAgenelist2           <- as.vector(CD20_elo_GSEAgenelist$stBeta)
names(CD20_elo_GSEAgenelist2)    <- CD20_elo_GSEAgenelist$geneSymbol

## run and plot fgsea for each cell type for agonism asymmetry and elo
fgseaRes_CD4_agAsym <- fgsea(pathways = pathways.hallmark,
                             stats = CD4_agAsym_GSEAgenelist2,
                             eps = 0.0)

fgseaRes_CD4_elo <- fgsea(pathways = pathways.hallmark,
                          stats = CD4_elo_GSEAgenelist2,
                          eps = 0.0)

fgseaRes_CD8_agAsym <- fgsea(pathways = pathways.hallmark,
                             stats = CD8_agAsym_GSEAgenelist2,
                             eps = 0.0)

fgseaRes_CD8_elo <- fgsea(pathways = pathways.hallmark,
                          stats = CD8_elo_GSEAgenelist2,
                          eps = 0.0)

fgseaRes_CD16_agAsym <- fgsea(pathways = pathways.hallmark,
                              stats = CD16_agAsym_GSEAgenelist2,
                              eps = 0.0)

fgseaRes_CD16_elo <- fgsea(pathways = pathways.hallmark,
                           stats = CD16_elo_GSEAgenelist2,
                           eps = 0.0)

fgseaRes_CD20_agAsym <- fgsea(pathways = pathways.hallmark,
                              stats = CD20_agAsym_GSEAgenelist2,
                              eps = 0.0)

fgseaRes_CD20_elo <- fgsea(pathways = pathways.hallmark,
                           stats = CD20_elo_GSEAgenelist2,
                           eps = 0.0,
                           minSize  = 0)
