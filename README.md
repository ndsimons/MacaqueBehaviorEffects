<p align="center">
  <img src="./misc/bx0bwn.jpg" width="700" title="hover text">
</p>
photo credit: TAO Photography

# Agonism and grooming behavior explain social status effects on physiology and gene regulation in rhesus macaques // Simons et al. 2021

This repository contains data and code for the manuscript _Agonism and grooming behavior explain social status effects on physiology and gene regulation in rhesus macaques - in prep_

#### All scripts can be run independently, but the following order follows the order of the manuscript:
##### 1. hierarchyMetrics.R - code to calculate dominance hierarchy metrics
##### 2. flowSortedModeling.R - model behavioral and elo effects in flow-sorted gene expression data set
##### 3. LPS_NULLmodeling.R - model behavioral and elo effects in LPS-challenged and NULL gene expression data
##### 4. gardModeling.R - model behavioral and elo effects in GARD-challenged gene expression data
##### 5. DEXmodeling.R - model behavioral and elo effects in DEX-challenged gene expression and chromatin accessibility data
##### 6. glucocorticoid_modeling.R - model behavioral and elo effects on glucocorticoid phenotypes (serum cortisol)
##### 7. mtDNA_modeling.R - model behavioral and elo effects on mtDNA copy number
##### 8. fgseaCode.R - code to run fgsea on elo & agonism asymmetry flow-sorted model results
##### 9. agAsymModel.R - model the effects of agonism asymmetry (residuals) in residual gene expression variance not explained by elo

