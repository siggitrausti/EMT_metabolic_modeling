# Script to score the essential genes for the D492M (MES) metabolism based on proteomic data
# Here, we use the METABRIC breast cancer data (https://www.cbioportal.org/study/summary?id=brca_metabric) 
# to identify the genes that are the best biomarkers for patients with claudin-low breast cancer, as these
# have been shown to display a EMT gene expression signature

# Sigurdur Karvelsson

library(siggitRausti)
library(survminer)
library(survival)
library(survcomp)
library(rmeta)
library(forestplot)

load('Dataset_2_130519.Rda') # DESeq2 normalized counts from the data
load('Dataset_2_clinical_130519.Rda') # Clinical information
clin_data <- clin_dat2
claudin_low_patients <- clin_dat2$PATIENT_ID[which(clin_dat2$CLAUDIN_SUBTYPE == 'claudin-low')]

# Essential genes:
MES_essential_genes <- c('ASL','IDH2','GUK1','PDHX','PRODH','PGLS','OAT','RENBP','CTH')

# Create a dataframe comprised of survival time, event and expression of the genes of interest:
surv_data <- clin_data[which(clin_data$PATIENT_ID %in% claudin_low_patients),]
exp_mat <- data.frame(t(MMS_rna_data3[,-c(1,2)]))
idx_vec <- rep(NA,nrow(surv_data))
for (i in 1:nrow(surv_data)){
  temp_id <- which(rownames(exp_mat) %in% surv_data$PATIENT_ID[i])
  if (length(temp_id) != 0){
    idx_vec[i] <- temp_id
  }
}

surv_data$ASL <- exp_mat$ASL[idx_vec]
surv_data$OAT <- exp_mat$OAT[idx_vec]
surv_data$PDHX <- exp_mat$PDHX[idx_vec]
surv_data$PRODH <- exp_mat$PRODH[idx_vec]
surv_data$RENBP <- exp_mat$RENBP[idx_vec]
surv_data$SLC19A3 <- exp_mat$SLC19A3[idx_vec]
surv_data$SLC7A5 <- exp_mat$SLC7A5[idx_vec]
surv_data$GUK1 <- exp_mat$GUK1[idx_vec]
surv_data$CTH <- exp_mat$CTH[idx_vec]
surv_data$IDH2 <- exp_mat$IDH2[idx_vec]
surv_data$PGLS <- exp_mat$PGLS[idx_vec]


surv_data$status <- ifelse(surv_data$OS_STATUS == 'LIVING',0,1)
surv_data$survival <- surv_data$OS_MONTHS

res.cut <- surv_cutpoint(surv_data, time = "survival", event = "status",
                         variables = MES_essential_genes)

summary(res.cut)
res.cat <- surv_categorize(res.cut)
res.cat <- res.cat[complete.cases(res.cat),]

# Prepare for the forest plot:
r.lower = rep(NA,length(MES_essential_genes)+1)
r.upper = rep(NA,length(MES_essential_genes)+1)
r.mean = rep(NA,length(MES_essential_genes)+1)
pvalues = rep(NA,length(MES_essential_genes))


for (i in 1:length(MES_essential_genes)){
  res.cat2 <- res.cat[complete.cases(res.cat),]
  y <- 'Surv(survival, status)'
  x <- MES_essential_genes[i]
  form = as.formula(paste(y, "~", x))
  # Generate a univariate cox proportional hazard model to use for C-index calculations:
  coxph_temp <- coxph(formula = form , data = res.cat)
  
  
  test_pred <- predict(coxph_temp,newdata = res.cat2, type = "lp")
  cindex_validation = concordance.index (test_pred, surv.time = res.cat2$survival,
                                         surv.event=res.cat2$status, method = "noether",alpha = 0.05)
  
  r.lower[i+1] <- cindex_validation$lower
  r.upper[i+1] <- cindex_validation$upper
  r.mean[i+1] <- cindex_validation$c.index
  pvalues[i] <- cindex_validation$p.value
  rm(test_pred)
  rm(cindex_validation)
  rm(coxph_temp)
  #rm(data_temp)
}
pvals_adj <- pvalues*length(pvalues)
pvals_adj <- p.adjust(pvalues,method = 'bonferroni')

#png('HER2_signature_forestplot_JAN2020.png',res=300,width = 6,height = 2.1,units = 'in')
tabletext <- cbind(c("Gene", "ASL", "IDH2", 
                     "GUK1","PDHX","PRODH","PGLS","OAT","RENBP","CTH"),c('Adj. p-value',paste0(signif(pvals_adj,2),c(' *',' ',' ',' ',' ',' ',' ',' ',' '))))

#png('Concordance_genes_DEC2020.png',res=300,width = 7,height = 4,units = 'in')
forestplot(tabletext,graph.pos = 2,lty.ci = 2,line.margin = .1,zero = 0.5,is.summary=c(TRUE,rep(FALSE,9)),
           boxsize = .25,mean=r.mean, lower=r.lower, upper=r.upper, 
           clip =c(0, 1),xticks = seq(0.5,1.0,0.1), 
           xlab="Concordance Index",col = fpColors(box=c("darkred")),
           title = 'C-index of essential proteomic MES genes')
#dev.off()