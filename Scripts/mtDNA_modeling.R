###################################################################
#### Code for modeling behavioral effects on mtDNA copy number ####
###################################################################
library(lmerTest)

## read in metadata
mtDNA_metadata <- read.delim('./Data/mtDNA_metadata')

## create results matrix
mtDNA_modelResults <- matrix(nrow = 13, ncol = 5)

## define list of behavioral variables for modeling 
behaviorList <- c('scaledElo','PC1.agg','PC1.groom','gc.noncontact.agg.rec','gc.contact.agg.rec','gc.noncontact.agg.given','gc.contact.agg.given','gc.groom.given','gc.groom.rec','gc.overall.groom','gc.overall.agg.rec','gc.overall.agg.given','gc.agg.Diff')
rownames(mtDNA_modelResults) <- behaviorList
colnames(mtDNA_modelResults) <- c('beta','SE(beta)','t','p','AIC')

## run linear mixed effects models and register results
for (i in behaviorList){
  tmpBehavior <- mtDNA_metadata[,i]
  mtDNA_modelResults[i,1:4] <- summary(lmer(logCN ~ tmpBehavior + cell + age + (1|batch), data = mtDNA_metadata, na.action = na.omit))$coefficients[c(2,9,23,30)]
}

## add AIC to results matrix
for (i in behaviorList){
  tmpBehavior <- mtDNA_metadata[,i]
  mtDNA_modelResults[i,5] <- AIC(lmer(logCN ~ tmpBehavior + cell + age + (1|batch), data = mtDNA_metadata, na.action = na.omit))
}

# get null AIC
AIC(lmer(logCN ~ cell + age + (1|batch), data = mtDNA_metadata, na.action = na.omit))
