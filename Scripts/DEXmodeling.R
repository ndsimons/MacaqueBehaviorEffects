##################################################################
#### DEX-challenge data set modeling (DEX GE, DEXΔ GE, DEX CA ####
##################################################################

## load libraries

library(limma)
library(edgeR)
library(EMMREML)
library(gridExtra)
library(grid)
library(cobs)
library(scales)
library(parallel)
library(doParallel)
library(ggplot2)
library(qvalue)

## load in data
load('./Data/dexData.RData')


#### Chromatin accessibility (DEX CA) model ####

## normalize the count matrix based on the number of reads mapped to the nuclear genome
voom_CA <- voom(calcNormFactors(DGEList(counts=ATAC_counts[,6:91],lib.size = ca_info$nuc_reads_q10)),plot=FALSE)

## Use limma to remove technical effects associated with the group (i.e., cage) and library prepartion (proportion of mtDNA-mapped reads)
design <- model.matrix(~ca_info$p_mtDNA+ca_info$group)
fit <-lmFit(voom_CA,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_CA=apply(residuals.MArrayLM(object=fit, voom_CA),2,function(x){x+intercepts})
rm(design); rm(fit)

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
model_ca <- function(behaveVar){
  design=model.matrix(~ca_info$trt+ca_info$TC1+ca_info$TC2+ca_info$TC3+ca_info$trt:ca_info[,behaveVar])
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix)
  res_full=resid_CA[,1:(3*ncol(design)+1)]
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
  colnames(res_full)[ncol(res_full)] = 'AIC'
  
  for(i in 1:nrow(resid_CA))
  {
    emma=emmreml(y=resid_CA[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    aic <- 2*5-2*emma$loglik
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"],aic))
  }
  res_full <- as.data.frame(res_full)
  #res_full$qval_NC <- qvalue(res_full[,20])$qvalue
  #res_full$qval_DEX <- qvalue(res_full[,21])$qvalue
  return(res_full)
}

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- model_ca(i)
  assign(value = behaveVar_tmp, x = paste(i,'_ca_modelRes_AIC',sep=''))
}

#### DEX CA permutations ####
## highly recommend running prmutations on a compute cluster

# Declare flag variable (1 at the first entry per animal, 0 in the rest, instrumental for reshuffling)
ca_info$flag <- rep(c(1,0), times=43)
ge_info$flag <- rep(c(1,0), times=43)

iters = 100
ca.perm <- function(behaveVar){
  ca_random<-ca_info
  for(iter in 1:iters)
  {
    print(iter)
    
    #Permute Elo among the set of samples flagged as 1 (one sample per individual)
    ca_random[,behaveVar][which(ca_random$flag==1)]=sample(ca_random[,behaveVar][which(ca_random$flag==1)])
    
    #Complete the permuted Elos of the other samples (this way the mapping individual-to-rank is conserved in the permutations)
    for(i in 1:length(levels(ca_random$shortID)))
    {
      set=which(ca_random$shortID==levels(ca_random$shortID)[i] & ca_random$flag==0)
      set_ref=which(ca_random$shortID==levels(ca_random$shortID)[i] & ca_random$flag==1)
      if(length(set)==1)
        ca_random[,behaveVar][set]=ca_random[,behaveVar][set_ref]
    }
    
    #Declare null (i.e. based on permutations) nested design for fixed effects.
    design = model.matrix(~ca_random$trt+ca_random$TC1+ca_random$TC2+ca_random$TC3+ca_random$trt:ca_random[,behaveVar])
    
    #Declare object res_null to store in it the permutations' p-values:
    res_null=resid_CA[,1:(ncol(design))]
    colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
    
    #Fit a model for each gene using emmreml
    for(i in 1:nrow(resid_CA))
    {
      emma=emmreml(y=resid_CA[i,],X=design,Z=as.matrix(Z_matrix),K=as.matrix(kin),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
      res_null[i,]=t(c(emma$pvalbeta[,"none"]))
    }
    
    #we register p-values of the associations to Elo at NC and LPS alone.
    if(iter==1)
    {
      shuffled_pvals_NC <-data.frame(x=res_null[,"p_value_ca_random$trtcontrol:ca_random[, behaveVar]"])
      shuffled_pvals_DEX <-data.frame(x=res_null[,"p_value_ca_random$trtdexamethasone:ca_random[, behaveVar]"])
      
      rownames(shuffled_pvals_NC)=rownames(res_null)
      rownames(shuffled_pvals_DEX)=rownames(res_null)
    } else {
      shuffled_pvals_NC <- cbind(shuffled_pvals_NC,x=res_null[,"p_value_ca_random$trtcontrol:ca_random[, behaveVar]"])
      shuffled_pvals_DEX <- cbind(shuffled_pvals_DEX,x=res_null[,"p_value_ca_random$trtdexamethasone:ca_random[, behaveVar]"])
    }
  }
}


for (i in behaviorList){
  print(i)
  behaveVar_tmp <- ca.perm(i)
  assign(value = behaveVar_tmp, x = paste(i,'_ca_perms',sep=''))
}


#### Gene expression (DEX GE) model ####
## voom normalize counts
voom_RNA <- voom(calcNormFactors(DGEList(counts=RNA_counts,lib.size = ge_info$mapped_reads)),plot=FALSE)

## remove technical effects associated with social group 
design <- model.matrix(~ge_info$group)
fit <-lmFit(voom_RNA,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_RNA=apply(residuals.MArrayLM(object=fit, voom_RNA),2,function(x){x+intercepts})
rm(design); rm(fit); rm(intercepts)

model_ge <- function(behaveVar){
  design=model.matrix(~ge_info$trt+ge_info$TC1+ge_info$TC2+ge_info$TC3+ge_info$trt:ge_info[,behaveVar])
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix)
  res_full=resid_RNA[,1:(3*ncol(design)+1)]
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
  colnames(res_full)[ncol(res_full)] = 'AIC'
  
  for(i in 1:nrow(resid_RNA))
  {
    emma=emmreml(y=resid_RNA[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    aic <- 2*5-2*emma$loglik
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"],aic))
  }
  res_full <- as.data.frame(res_full)
  return(res_full)
}

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- model_ge(i)
  assign(value = behaveVar_tmp, x = paste(i,'_ge_modelRes_AIC',sep=''))
}

#### DEX GE permutations ####
## highly recommend running prmutations on a compute cluster

# Declare flag variable (1 at the first entry per animal, 0 in the rest, instrumental for reshuffling)
ca_info$flag <- rep(c(1,0), times=43)
ge_info$flag <- rep(c(1,0), times=43)

iters = 100
ge.perm <- function(behaveVar){
  ge_random<-ge_info
  for(iter in 1:iters)
  {
    print(iter)
    
    #Permute Elo among the set of samples flagged as 1 (one sample per individual)
    ge_random[,behaveVar][which(ge_random$flag==1)]=sample(ge_random[,behaveVar][which(ge_random$flag==1)])
    
    #Complete the permuted Elos of the other samples (this way the mapping individual-to-rank is conserved in the permutations)
    for(i in 1:length(levels(ge_random$shortID)))
    {
      set=which(ge_random$shortID==levels(ge_random$shortID)[i] & ge_random$flag==0)
      set_ref=which(ge_random$shortID==levels(ge_random$shortID)[i] & ge_random$flag==1)
      if(length(set)==1)
        ge_random[,behaveVar][set]=ge_random[,behaveVar][set_ref]
    }
    
    #Declare null (i.e. based on permutations) nested design for fixed effects.
    design = model.matrix(~ge_random$trt+ge_random$TC1+ge_random$TC2+ge_random$TC3+ge_random$trt:ge_random[,behaveVar])
    
    #Declare object res_null to store in it the permutations' p-values:
    res_null=resid_RNA[,1:(ncol(design))]
    colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
    
    #Fit a model for each gene using emmreml
    for(i in 1:nrow(resid_RNA))
    {
      emma=emmreml(y=resid_RNA[i,],X=design,Z=as.matrix(Z_matrix),K=as.matrix(kin),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
      res_null[i,]=t(c(emma$pvalbeta[,"none"]))
    }
    
    #we register p-values of the associations to Elo at NC and LPS alone.
    if(iter==1)
    {
      shuffled_pvals_NC <-data.frame(x=res_null[,"p_value_ge_random$trtcontrol:ge_random[, behaveVar]"])
      shuffled_pvals_DEX <-data.frame(x=res_null[,"p_value_ge_random$trtdex:ge_random[, behaveVar]"])
      
      rownames(shuffled_pvals_NC)=rownames(res_null)
      rownames(shuffled_pvals_DEX)=rownames(res_null)
    } else {
      shuffled_pvals_NC <- cbind(shuffled_pvals_NC,x=res_null[,"p_value_ge_random$trtcontrol:ge_random[, behaveVar]"])
      shuffled_pvals_DEX <- cbind(shuffled_pvals_DEX,x=res_null[,"p_value_ge_random$trtdex:ge_random[, behaveVar]"])
    }
  }
}


for (i in behaviorList){
  print(i)
  behaveVar_tmp <- ge.perm(i)
  assign(value = behaveVar_tmp, x = paste(i,'_ge_perms',sep=''))
}


#### Δ Gene expression (DEX GEΔ) model ####
voom_RNA_delta=voom_RNA$E[,seq(1,86,2)]-voom_RNA$E[,seq(2,86,2)]
design <- model.matrix(~ge_info[seq(1,86,2),"group"])
fit <-lmFit(voom_RNA_delta,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_RNA_delta=apply(residuals.MArrayLM(object=fit, voom_RNA_delta),2,function(x){x+intercepts})
rm(design); rm(fit); rm(intercepts)

ge_info_delta <- ge_info[seq(1,86,2),-2]

Z_matrix_delta <- Z_matrix[seq(from=1,to=86,by=2),]

model_geDelta <- function(behaveVar){
  design=model.matrix(~ge_info_delta$TC1+ge_info_delta$TC2+ge_info_delta$TC3+ge_info_delta[,behaveVar])
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix[seq(1,86,2),])
  res_full=resid_RNA_delta[,1:(3*ncol(design)+1)]
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
  colnames(res_full)[ncol(res_full)] = 'AIC'
  
  for(i in 1:nrow(resid_RNA))
  {
    emma=emmreml(y=resid_RNA_delta[i,],X=design,Z=as.matrix(Z_matrix_delta),K=as.matrix(kin),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    aic <- 2*4-2*emma$loglik
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"],aic))
  }
  res_full <- as.data.frame(res_full)
  return(res_full)
}

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- model_geDelta(i)
  assign(value = behaveVar_tmp, x = paste(i,'_GE_delta_modelRes_AIC',sep=''))
}

#### DEX GE delta  permutations ####
## highly recommend running prmutations on a compute cluster


iters = 2
ge_delta.perm <- function(behaveVar){
  ge_random<-ge_info_delta
  for(iter in 1:iters)
  {
    print(iter)
    
    #Declare null (i.e. based on permutations) nested design for fixed effects.
    randVar <- sample(ge_random[,behaveVar])
    design = model.matrix(~ge_random$TC1+ge_random$TC2+ge_random$TC3+randVar)
    #Declare object res_null to store in it the permutations' p-values:
    res_null=resid_RNA_delta[,1:(ncol(design))]
    colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
    
    #Fit a model for each gene using emmreml
    for(i in 1:nrow(resid_RNA_delta))
    {
      emma=emmreml(y=resid_RNA_delta[i,],X=design,Z=as.matrix(Z_matrix_delta),K=as.matrix(kin),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
      res_null[i,]=t(c(emma$pvalbeta[,"none"]))
    }
    #we register p-values of the associations to Elo at NC and LPS alone.
    if(iter==1)
    {
      shuffled_pvals_ge_delta <-data.frame(x=res_null[,"p_value_randVar"])
      rownames(shuffled_pvals_ge_delta)=rownames(res_null)
    } else {
      shuffled_pvals_ge_delta <- cbind(shuffled_pvals_ge_delta,x=res_null[,"p_value_randVar"])
    }
  }}
  

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- ge.perm(i)
  assign(value = behaveVar_tmp, x = paste(i,'_ge_delta_perms',sep=''))
}
