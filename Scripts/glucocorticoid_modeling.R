###########################################################################
#### Code for modeling behavioral effects on glucocorticoid phenotypes ####
###########################################################################

library(lmerTest)

## read in cortisol data and metadata
cort_data <- read.table('/.Data/CortData.txt', header =T)
cort_metadata <- read.table('/.Data/cort_metadata.txt', header =T)

## assign list of behavioral variables
behaviorList <- c('scaledElo',
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
                  "gc.agg.Diff",
                  'NULL')


#### Diurnal serum cortisol modeling ####
cort_data_diurnal <- reshape2::melt(cort_data[,1:3])
cort_metadata_diurnal <- subset(cort_metadata, Phase == 1)

gc_diurnal_resMatrix <- matrix(nrow = 14, ncol = 5)
rownames(gc_diurnal_resMatrix)[1:14] <- behaviorList
colnames(gc_diurnal_resMatrix) <- c('beta','SE(beta)','t','p','AIC')

## run linear mixed effects models and register results
for (i in behaviorList){
  tmpBehavior <- cort_metadata_diurnal[,i]
  gc_diurnal_resMatrix[i,1:4] <- summary(lmer(cort_data_diurnal$value ~ tmpBehavior + cort_metadata_diurnal$Time + cort_metadata_diurnal$Age + (1|cort_metadata_diurnal$ID)))$coefficients[c(2,6,14,18)]
}

## add AIC to results matrix
for (i in behaviorList){
  tmpBehavior <- cort_metadata_diurnal[,i]
  gc_diurnal_resMatrix[i,5] <- AIC(lmer(cort_data_diurnal$value ~ tmpBehavior + cort_metadata_diurnal$Time + cort_metadata_diurnal$Age + (1|cort_metadata_diurnal$ID)))
}

## get NULL model AIC
gc_diurnal_resMatrix[14,5] <- AIC(lmer(cort_data_diurnal$value ~ cort_metadata_diurnal$Time + cort_metadata_diurnal$Age + (1|cort_metadata_diurnal$ID)))


#### cortisol slope modeling ####
cort_slope <- cort_data[,7]
cort_dct_1.5 <- cort_data[,10]
cort_dct_4.5 <- cort_data[,13]
cort_dst <- cort_data[,16]

cort_metadata_sub <- subset(cort_metadata, Time == '8' & Phase == '1')
cort_metadata_sub$control_cort <- c(cort_data$X800_P1)
cort_metadata_sub$control_dct <- c(cort_data$p1_dct_baseline)
cort_metadata_sub$control_dst <- c(cort_data$X1100_P1)
cort_metadata_sub$control_dst_dex <- c(cort_data$p1_dex_dst)
cort_metadata_sub$control_dct_dex_1.5 <- c(cort_data$p1_dex_1.5)
cort_metadata_sub$control_dct_dex_4.5 <- c(cort_data$p1_dex_4.5)

cort_slope_resMatrix <- matrix(nrow = 14, ncol = 5)
rownames(cort_slope_resMatrix)[1:14] <- behaviorList
colnames(cort_slope_resMatrix) <- c('beta','SE(beta)','t','p','AIC')
## run linear mixed effects models and register results
for (i in behaviorList){
  tmpBehavior <- cort_metadata_sub[,i]
  cort_slope_resMatrix[i,1:4] <- summary(lm(cort_slope ~ tmpBehavior + cort_metadata_sub$control_cort + cort_metadata_sub$Age))$coefficients[c(2,6,10,14)]
}

## add AIC to results matrix
for (i in behaviorList){
  tmpBehavior <- cort_metadata_sub[,i]
  cort_slope_resMatrix[i,5] <- AIC(lm(cort_slope ~ tmpBehavior + cort_metadata_sub$control_cort + cort_metadata_sub$Age))
}

## get NULL model AIC
cort_slope_resMatrix[14,5] <- AIC(lm(cort_slope ~ cort_metadata_sub$control_cort + cort_metadata_sub$Age))


#### cort 1.5hr modeling ####
cort_dct_1.5_resMatrix <- matrix(nrow = 14, ncol = 5)
rownames(cort_dct_1.5_resMatrix)[1:14] <- behaviorList
colnames(cort_dct_1.5_resMatrix) <- c('beta','SE(beta)','t','p','AIC')
## run linear mixed effects models and register results
for (i in behaviorList){
  tmpBehavior <- cort_metadata_sub[,i]
  cort_dct_1.5_resMatrix[i,1:4] <- summary(lm(cort_dct_1.5 ~ tmpBehavior + tmp_metadata$control_dct + tmp_metadata$control_dct_dex_1.5 + tmp_metadata$Age))$coefficients[c(2,7,12,17)]
}

## add AIC to results matrix
for (i in behaviorList){
  tmpBehavior <- cort_metadata_sub[,i]
  cort_dct_1.5_resMatrix[i,5] <- AIC(lm(cort_dct_1.5 ~ tmpBehavior + tmp_metadata$control_dct + tmp_metadata$control_dct_dex_1.5 + tmp_metadata$Age))
}

## get NULL model AIC
cort_dct_1.5_resMatrix[14,5] <- AIC(lm(cort_dct_1.5 ~ tmp_metadata$control_dct + tmp_metadata$control_dct_dex_1.5 + tmp_metadata$Age))


#### cort 4.5hr modeling ####
cort_dct_4.5_resMatrix <- matrix(nrow = 14, ncol = 5)
rownames(cort_dct_4.5_resMatrix)[1:14] <- behaviorList
colnames(cort_dct_4.5_resMatrix) <- c('beta','SE(beta)','t','p','AIC')
## run linear mixed effects models and register results
for (i in behaviorList){
  tmpBehavior <- cort_metadata_sub[,i]
  cort_dct_4.5_resMatrix[i,1:4] <- summary(lm(cort_dct_4.5 ~ tmpBehavior + tmp_metadata$control_dct + tmp_metadata$control_dct_dex_4.5 + tmp_metadata$Age))$coefficients[c(2,7,12,17)]
}

## add AIC to results matrix
for (i in behaviorList){
  tmpBehavior <- cort_metadata_sub[,i]
  cort_dct_4.5_resMatrix[i,5] <- AIC(lm(cort_dct_4.5 ~ tmpBehavior + tmp_metadata$control_dct + tmp_metadata$control_dct_dex_4.5 + tmp_metadata$Age))
}

## get NULL model AIC
cort_dct_4.5_resMatrix[14,5] <- AIC(lm(cort_dct_4.5 ~ tmp_metadata$control_dct + tmp_metadata$control_dct_dex_4.5 + tmp_metadata$Age))


#### cort 24hr modeling ####
cort_dst_resMatrix <- matrix(nrow = 14, ncol = 5)
rownames(cort_dst_resMatrix)[1:14] <- behaviorList
colnames(cort_dst_resMatrix) <- c('beta','SE(beta)','t','p','AIC')
## run linear mixed effects models and register results
for (i in behaviorList){
  tmpBehavior <- cort_metadata_sub[,i]
  cort_dst_resMatrix[i,1:4] <- summary(lm(cort_dst ~ tmpBehavior + tmp_metadata$control_dst + tmp_metadata$control_dst_dex + tmp_metadata$Age ))$coefficients[c(2,7,12,17)]
}

## add AIC to results matrix
for (i in behaviorList){
  tmpBehavior <- cort_metadata_sub[,i]
  cort_dst_resMatrix[i,5] <- AIC(lm(cort_dst ~ tmpBehavior + tmp_metadata$control_dst + tmp_metadata$control_dst_dex + tmp_metadata$Age ))
}

## get NULL model AIC
cort_dst_resMatrix[14,5] <- AIC(lm(cort_dct_4.5 ~ tmp_metadata$control_dct + tmp_metadata$control_dct_dex_4.5 + tmp_metadata$Age))
