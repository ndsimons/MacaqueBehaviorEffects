#############################################
#### GARD gene expression modeling ####
#############################################

## load libraries
library(limma)
library(edgeR)
library(EMMREML)
library(nlme)
library(qvalue)
library(openxlsx)
library(cobs)


## load GARD data, metadata, kinship matrix, and Z matrix
load('./Data/gardData.RData')

## Use limma to remove technical effects associated with the group 
design <- model.matrix(~Gard_metadata$Group)
voom_Gard <- voom(calcNormFactors(DGEList(counts=Gard_filtered_reads)),plot=FALSE)
fit <-lmFit(voom_Gard,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_RNA=apply(residuals.MArrayLM(object=fit, voom_Gard),2,function(x){x+intercepts})
rm(design); rm(fit); rm(intercepts)

#### Run gard models ####
behaviorList <- c('Elo',
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
model_gard <- function(behaveVar){
  design=model.matrix(~Gard_metadata[,behaveVar]+Gard_metadata$Age+Gard_metadata$PC1.TC+Gard_metadata$PC2.TC)
  K=as.matrix(K)
  Z=as.matrix(Z)
  res_full=resid_RNA[,1:(3*ncol(design)+1)]
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
  colnames(res_full)[ncol(res_full)] = 'AIC'
  
  for(i in 1:nrow(resid_RNA))
  {
    emma=emmreml(y=resid_RNA[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    aic <- 2*4-2*emma$loglik
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"],aic))
  }
  res_full <- as.data.frame(res_full)
  #res_full$qval_NC <- qvalue(res_full[,20])$qvalue
  #res_full$qval_DEX <- qvalue(res_full[,21])$qvalue
  return(res_full)
}

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- model_gard(i)
  assign(value = behaveVar_tmp, x = paste(i,'_gard_modelRes',sep=''))
}


#### GARD permutations ####
## highly recommend running permutations on a compute cluster

iters = 100

gard.perm <- function(behaveVar){
  cols_random<-Gard_metadata
  for(iter in 1:iters)
  {
    print(iter)

    #Declare design
    design=model.matrix(~sample(cols_random[,behaveVar])+cols_random$Age+cols_random$PC1.TC+cols_random$PC2.TC)
    
    #Declare object res_null to store in it the permutations' p-values:
    res_null=resid_RNA[,1:(ncol(design))]
    colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
    
    #Fit a model for each gene using emmreml
    for(i in 1:nrow(resid_RNA))
    {
      emma=emmreml(y=resid_RNA[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
      res_null[i,]=t(c(emma$pvalbeta[,"none"]))
    }
    
    #we register p-values of the associations to Elo at NC and LPS alone.
    if(iter==1)
    {
      shuffled_pvals <-data.frame(x=res_null[,"p_value_sample(cols_random[, behaveVar])"])
      rownames(shuffled_pvals)=rownames(res_null)
    } else {
      shuffled_pvals <- cbind(shuffled_pvals,x=res_null[,"p_value_sample(cols_random[, behaveVar])"])
    }
  }
  return(shuffled_pvals)
}

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- gard.perm(i)
  assign(value = behaveVar_tmp, x = paste(i,'_gard_perms',sep=''))
}

