##############################################
#### Flow-sorted gene expression modeling ####
##############################################

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
load('./Data/flowSortedData.RData')

## Voom normalize and remove group ('batch') effects using limma:
residuals=lapply(names(info_cell),function(x){
  design <- model.matrix(~0+info_cell[[x]]$group)
  dge <- DGEList(counts=read_count_per_cell[[x]])
  dge <- calcNormFactors(dge)
  v <- voom(dge,design,plot=FALSE)
  fit <-lmFit(v,design)
  return(residuals.MArrayLM(object=fit, v))
});names(residuals)=names(info_cell)

#### Models ####
cellsort.Model <- function(behaveVar){
  lapply(c("CD4","CD8","CD14","CD16","CD20"),function(x){
    tmp=t(apply(residuals[[x]],1,function(y){
      design=model.matrix(~info_cell[[x]][,behaveVar]+info_cell[[x]]$age)
      K=as.matrix(kin_per_cell[[x]])
      Z=as.matrix(Z_matrix[[x]])
      emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
      p=emma$pvalbeta ## store p-values
      varb=emma$varbetahat ## store variance of BLUPs
      b=emma$betahat ## store BLUPs
      aic <- 2*3-2*emma$loglik # calculate AIC for each model
      return(c(b,varb,p[,"none"],aic)) ## return results for each gene
    })); colnames(tmp)=c('beta_intercept','beta_behaveVar','beta_age','var_intercept','var_behaveVar','var_age','pval_intercept','pval_behaveVar','pval_age','AIC'); tmp <- as.data.frame(tmp); tmp$qval <- qvalue(tmp$pval_behaveVar)$qvalue;return(tmp)})
}

behaviorList <- c('elo',
                  'PC1.agg',
                  'PC1.groom',
                  "gc.noncontact.agg.rec",
                  "gc.contact.agg.rec",
                  "gc.noncontact.agg.given",
                  "gc.contact.agg.given",
                  "gc.groom.given",
                  "gc.groom.rec",
                  "gc.overall.groom",
                  "gc.overall.agg.rec",
                  "gc.overall.agg.given",
                  "gc.agg.Diff")

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- cellsort.Model(i)
  names(behaveVar_tmp) <- c("CD4","CD8","CD14","CD16","CD20")
  saveRDS(behaveVar_tmp, file=paste0(i,'_modelResAIC.rds'))
  assign(value = behaveVar_tmp, x = paste(i,'_modelResAIC',sep=''))
}

######################
#### Permutations ####
######################

## highly recommended to run this on a computing cluster
perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
  
  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df[order(input_df[,pvals_index]),]
  p_obs <- data.frame(pvalue=ro[,pvals_index])
  p_vector<-matrix(as.matrix(perm_df),ncol=1)
  p_vector=data.frame(p_vector[order(p_vector)])
  
  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]
  
  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])
  
  for(i in 1:observed)
  {
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
    if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
  }
  tabla <-data.frame(pi_hat,pval=p_obs[,1])
  
  tabla[1,]=c(1,0)
  last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
  tabla[nrow(tabla),]=c(last_percentile_average,1)
  constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
  f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))
  
  f_hat_serie=f_hat$fitted
  pi_o=f_hat_serie[length(f_hat_serie)]
  pi_o=min(pi_o,1)
  pi_o=max(pi_o,0)
  
  Fdr_ST_perm=pi_o*F_o/F
  
  for(i in 1:length(p_obs[,1]))
  {
    Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
      }
    }
    if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
  }
  
  fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
  rownames(fdrs_df)=rownames(ro)
  colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)
  
  return(fdrs_df)
}
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

#interate through each behavior and calculate permuted p vals for each cell type
for (i in behaviorList){
  print(i)
  behaveVar_tmp <- cellsort.Perm(i)
  names(behaveVar_tmp) <- c("CD4","CD8","CD14","CD16","CD20")
  saveRDS(behaveVar_tmp, file=paste0(i,'_perms.rds'))
  assign(value = behaveVar_tmp, x = paste(i,'_perms',sep=''))
}


# iterate through each behavior and calculate the fdrs and return object with fdrs added in
for (i in behaviorList){
  print(i)
  tmp <- paste(i,'_modelResAIC',sep='')
  res <- get(tmp)
  tmp <- paste(i,'_perms',sep='')
  perms <- get(tmp)
  for (n in names(res)){
    res[[n]] <- perm.fdr(res[[n]],perms[[n]][,2],'pval_behaveVar','behaveVar')
  }
  saveRDS(res, file=paste0(i,'_FDR.rds'))
  assign(value=res, x=paste(i,'_modelResAIC_FDR',sep=''))
}
