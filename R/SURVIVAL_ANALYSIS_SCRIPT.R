# Script to perform the survival analysis on the metadata of the patients from where the proteomic data was acquired

# Sigurdur Karvelsson

library(siggitRausti)
library(tidyr)
library(survminer)
library(survival)
library(readxl)
library(dplyr)

# Now perform the survival analysis of the patients:
# read in data (Tang et al (2018).
tang <- read_excel('13073_2018_602_MOESM3_ESM (2).xlsx','PSM')
tang_meta <- read_excel('13073_2018_602_MOESM3_ESM (2).xlsx','patient_characteristics',col_names = F)
colnames(tang_meta)[c(1:3)] <- c('patient','subtype','ER')
tang_meta$death <- NA
tang_meta$survival <- NA
tang_meta$age <- NA
tang_meta$ER2 <- NA

# Get the survival information:
library(GEOquery)
gse <- getGEO("GSE37751", GSEMatrix = TRUE)
names(pData(gse[[1]]))

for (i in 1:nrow(tang_meta)){
  id_in_meta <- grep(paste('',substring(tang_meta$patient[i],4)),pData(gse[[1]])[,24])
  print(id_in_meta)
  if (length(id_in_meta) != 0){
    tang_meta$death[i] <- pData(gse[[1]])[,47][id_in_meta]
    tang_meta$survival[i] <- pData(gse[[1]])[,54][id_in_meta]
    tang_meta$age[i] <- pData(gse[[1]])[,45][id_in_meta]
    tang_meta$ER2[i] <- pData(gse[[1]])[,48][id_in_meta]
  }
}

tang_meta$death[is.na(tang_meta$death)] <- 'not applicable'
for (i in 1:nrow(tang_meta)){
  if (tang_meta$death[i] == 'not applicable'){
    tang_meta$death[i] <- NA
  } else if (tang_meta$death[i] == 'Yes'){
    tang_meta$death[i] <- 1
  } else if (tang_meta$death[i] == 'No'){
    tang_meta$death[i] <- 0
  } else if (is.na(tang_meta$death[i])){
    tang_meta$death[i] <- NA
  }
}
tang_meta$survival <- as.numeric(tang_meta$survival)
tang_meta$death <- as.numeric(tang_meta$death)
tang_meta$age <- as.numeric(tang_meta$age)

# Process the data a bit
tang_t <- data.frame(t(tang[,-c(1:3)]))
colnames(tang_t) <- tang$Symbol

# Introduce ASL proteomic levels:
tang_meta$ASL <- tang_t$ASL
genes_of_interest <- c('ASL')


# Create a dataset for ERpos and ERneg patients:
erpos_tang <- tang_meta[intersect(which(tang_meta$ER %in% 'ER+'),which(tang_meta$subtype %in% 'Tumor')),]
erneg_tang <- tang_meta[intersect(which(tang_meta$ER %in% 'ER-'),which(tang_meta$subtype %in% 'Tumor')),]


# Perform the univariate Cox-proportional hazard modeling for ERneg patients:
res.cox.neg.uni <- coxph(Surv(survival, death) ~ ASL, data =  erneg_tang)
summary(res.cox.neg.uni)

# Perform the multivariate Cox-proportional hazard modeling for ERneg patients:
res.cox.neg.multi <- coxph(Surv(survival, death) ~ age + ASL, data =  erneg_tang)
summary(res.cox.neg.multi)

# Perform the univariate Cox-proportional hazard modeling for ERpos patients:
res.cox.pos.uni <- coxph(Surv(survival, death) ~ ASL, data =  erpos_tang)
summary(res.cox.pos.uni)

# Perform the multivariate Cox-proportional hazard modeling for ERpos patients:
res.cox.pos.multi <- coxph(Surv(survival, death) ~ age + ASL, data =  erpos_tang)
summary(res.cox.pos.multi)




