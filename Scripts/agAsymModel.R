####################################################
#### Agonism asymmetry gene expression modeling ####
####################################################

## load libraries
library(ggplot2)
library(limma)
library(edgeR)
library(EMMREML)
library(data.table)
library(cobs)
library(scales)
library(qvalue)
library(GGally)
library(pheatmap)
library(parallel)
library(doParallel)


## load in flow-sorted gene expression data, metadata, Z matrix, and kinship matrix
load('./Data/agAsymResidualModelData.RData')

## Voom normalize and remove group and rank effects
residuals=lapply(names(info_cell),function(x){
  design <- model.matrix(~0+info_cell[[x]]$group + info_cell[[x]]$elo)
  dge <- DGEList(counts=read_count_per_cell[[x]])
  dge <- calcNormFactors(dge)
  v <- voom(dge,design,plot=FALSE)
  fit <-lmFit(v,design)
  return(residuals.MArrayLM(object=fit, v))
});names(residuals)=names(info_cell)

## generate residual agonism asymmetry from agAsym~elo 
info_cell=lapply(info_cell,function(x){
  tmp=x
  tmp$aggDiffResiduals <- resid(lm(tmp$gc.agg.Diff ~ tmp$elo))
  return(tmp)
})

#### Models ####
cellsort.residual.Model <- function(behaveVar){
  lapply(c("CD4","CD8","CD14","CD16","CD20"),function(x){
    tmp=t(apply(residuals[[x]],1,function(y){
      design=model.matrix(~info_cell[[x]][,behaveVar]+info_cell[[x]]$age)
      K=as.matrix(kin_per_cell[[x]])
      Z=as.matrix(Z_matrix[[x]])
      emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
      p=emma$pvalbeta ## store p-values
      varb=emma$varbetahat ## store variance of BLUPs
      b=emma$betahat ## store BLUPs
      aic <- 2*3-2*emma$loglik
      return(c(b,varb,p[,"none"],aic)) ## return results for each gene
    })); colnames(tmp)=c('beta_intercept','beta_behaveVar','beta_age','var_intercept','var_behaveVar','var_age','pval_intercept','pval_behaveVar','pval_age','AIC'); tmp <- as.data.frame(tmp); tmp$qval <- qvalue(tmp$pval_behaveVar)$qvalue;return(tmp)})
}


behaviorList <- c('aggDiffResiduals')
for (i in behaviorList){
  print(i)
  behaveVar_tmp <- cellsort.residual.Model(i)
  names(behaveVar_tmp) <- c("CD4","CD8","CD14","CD16","CD20")
  assign(value = behaveVar_tmp, x = paste(i,'_residualModelRes',sep=''))
}


#### Permutations ####
## highly recommend running permuations on a compute cluster
clus <- makeCluster(32)
registerDoParallel(cores=32)

cellsort.Perm <- function(behaveVar){
  EMMA_permute_rank_pvalues=lapply(names(residuals),function(pp){
    nperm=100
    for (sim in 1:nperm){
      cell=pp
      clusterExport(clus,c("residuals","info_cell","Z_matrix","kin_per_cell"))
      behaveVar_perm=sample(info_cell[[cell]][,behaveVar])
      tmp=parApply(clus,residuals[[cell]],1,function(y){
        .libPaths( c("/data/tunglab/nds/Rlibs/", .libPaths() ) )
        library(EMMREML)
        design=model.matrix(~behaveVar_perm+info_cell[[cell]]$age) 
        K=as.matrix(kin_per_cell[[cell]])
        Z=as.matrix(Z_matrix[[cell]])
        emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
        return(emma$pvalbeta[,"none"])}) 
      print(sim)
      if (sim==1) final=t(tmp)
      else final=rbind(final,t(tmp))
    }
    colnames(final)=c("intercept_pvalue","permuted_behaveVar_pvalue","age_pvalue")
    return(final)
  }
  )
}

## iterate through each behavior and calculate permuted p vals for each cell type
for (i in behaviorList){
  print(i)
  behaveVar_tmp <- cellsort.Perm(i)
  names(behaveVar_tmp) <- c("CD4","CD8","CD14","CD16","CD20")
  assign(value = behaveVar_tmp, x = paste(i,'_perms',sep=''))
}
