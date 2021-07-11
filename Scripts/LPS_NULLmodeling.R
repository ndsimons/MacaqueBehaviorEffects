#############################################
#### LPS & NULL gene expression modeling ####
#############################################

## load libraries
library(ggplot2)
library(limma)
library(edgeR)
library(EMMREML)
library(data.table)
library(cobs)

## Load input files: reads, metadata and kinship matrixes in three data.frames (reads, cols, K, respectively)
load('./Data/lpsData.R')

#### Run linear model: (nested in LPS+ and LPS-) ####

#Get dge object
dge <- DGEList(counts=reads)
dge <- calcNormFactors(dge)

#Remove group effects using limma and extract residuals.
design = model.matrix(~0+Group,data=cols)
v <- voom(dge,design,plot=FALSE)
fit <-lmFit(v,design)
fit <- eBayes(fit)
exp<-residuals.MArrayLM(object=fit, v)

#order genes alphabetically to ensure a criterium for ordering rows all the time.
exp=exp[order(rownames(exp)),]

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

model_lps <- function(behaveVar){
  design = model.matrix(~Condition+Condition:cols[,behaveVar]+Condition:Age+Condition:PC1+Condition:PC2,data=cols)
  
  #Declare object res_full to store in it the model results: beta coefficients, standard deviations and p values
  res_full=exp[,1:(3*ncol(design))]
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
  
  #Declare object random_effects to store in it the individual-wise u random effects
  random_effects=exp[,1:length(levels(cols$Animal))]
  colnames(random_effects)=levels(cols$Animal)
  
  #Define matrix Z describing the sample-to-individual mapping
  Z=matrix(rep(0,nrow(cols)*ncol(K)),nrow=nrow(cols),ncol=ncol(K))
  rownames(Z)=rownames(cols)
  colnames(Z)=colnames(K)
  for(i in 1:ncol(Z))
  {
    set=which(cols$Animal == colnames(Z)[i])
    Z[set,i]=1
  }
  
  #Fit a model for each gene using emmreml
  for(i in 1:nrow(exp))
  {
    emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    random_effects[i,]=t(emma$uhat)
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
  }
  res_full <- as.data.frame(res_full)
  return(res_full)
}

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- model_lps(i)
  assign(value = behaveVar_tmp, x = paste(i,'_lps_modelRes',sep=''))
}


#### LPS/NULL permutations ####
## highly recommend running permutations on a compute cluster

iters = 100

lps.perm <- function(behaveVar){
  cols_random<-cols
  for(iter in 1:iters)
  {
    print(iter)
    
    #Permute behaveVar among the set of samples flagged as 1 (one sample per individual)
    cols_random[which(cols_random$flag==1),behaveVar]=sample(cols_random[which(cols_random$flag==1),behaveVar])
    
    #Complete the permuted behaveVar of the other samples (this way the mapping individual-to-rank is conserved in the permutations)
    for(i in 1:length(levels(cols_random$Animal)))
    {
      set=which(cols_random$Animal==levels(cols_random$Animal)[i] & cols_random$flag==0)
      set_ref=which(cols_random$Animal==levels(cols_random$Animal)[i] & cols_random$flag==1)
      if(length(set)==1)
        cols_random[,behaveVar][set]=cols_random[,behaveVar][set_ref]
    }
    
    #Declare null (i.e. based on permutations) nested design for fixed effects.
    design = model.matrix(~Condition+Condition:cols_random[,behaveVar]+Condition:Age+Condition:PC1+Condition:PC2,data=cols_random)
    
    #Declare object res_null to store in it the permutations' p-values:
    res_null=exp[,1:(ncol(design))]
    colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
    
    #Fit a model for each gene using emmreml
    for(i in 1:nrow(exp))
    {
      emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
      res_null[i,]=t(c(emma$pvalbeta[,"none"]))
    }
    
    #we register p-values of the associations to Elo at NC and LPS alone.
    if(iter==1)
    {
      shuffled_pvals_NC <-data.frame(x=res_null[,"p_value_ConditionNC:cols_random[, behaveVar]"])
      shuffled_pvals_LPS <-data.frame(x=res_null[,"p_value_ConditionLPS:cols_random[, behaveVar]"])
      
      rownames(shuffled_pvals_NC)=rownames(res_null)
      rownames(shuffled_pvals_LPS)=rownames(res_null)
    } else {
      shuffled_pvals_NC <- cbind(shuffled_pvals_NC,x=res_null[,"p_value_ConditionNC:cols_random[, behaveVar]"])
      shuffled_pvals_LPS <- cbind(shuffled_pvals_LPS,x=res_null[,"p_value_ConditionLPS:cols_random[, behaveVar]"])
    }
  }
  return(shuffled_pvals_NC)
  return(shuffled_pvals_LPS)
}

for (i in behaviorList){
  print(i)
  behaveVar_tmp <- lps.perm(i)
  assign(value = behaveVar_tmp, x = paste(i,'_lps_perms',sep=''))
}


