# Script to perform the dataset correlation analysis

# Sigurdur Karvelsson


library(siggitRausti)
library(tidyr)
library(readxl)
library(dplyr)
library(biomaRt)
library(PerformanceAnalytics)

# Create a 3x3 correlation image with distributions. 

# Load all data:
rnaseq <- read_excel('D492_D492M_RNAseq.xlsx','ALL_USE_THIS')
prot <- read_excel('Analysis Results Qiong.xlsx')
microarray <- read_excel('E_vs_M.xls')
colnames(rnaseq) <- c('Genes','d492','d492m')
rnaseq$fold <- rnaseq$d492m/rnaseq$d492
rnaseq$type <- 'rnaseq'
rnaseq <- aggregate(fold ~ Genes, data = rnaseq, FUN = mean)
colnames(rnaseq)[2] <- 'D492M/D492'

# Work to get the proteomic data with gene names and D492M/D492 fold changes:
prot <- prot[,c(2,5,10,11,13,14,16)]
prot <- data.frame(prot)
prot$fold_change <- NA
for (i in 1:nrow(prot)){
  prot$fold_change[i] <- mean(as.numeric(prot[i,c(3,4,7)]))/mean(as.numeric(prot[i,c(2,5,6)]))
}
prot <- prot[,c(1,8)]
colnames(prot) <- c('Genes','temp')

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
hugo_genes <- getBM(attributes=c('hgnc_symbol','uniprot_gn_id','illumina_humanht_12_v3'), 
                    filters = 'uniprot_gn_id', 
                    #values = colnames(rnaseq_D492_2), 
                    values = prot$Genes,
                    mart = ensembl)
prot$Genes2 <- NA
for (i in 1:nrow(prot)){
  if (prot$Genes[i] %in% hugo_genes$uniprot_gn_id){
    prot$Genes2[i] <- unique(hugo_genes$hgnc_symbol[which(hugo_genes$uniprot_gn_id %in% prot$Genes[i])])
  }
}

prot <- prot[,c(3,2)]
prot <- prot[complete.cases(prot),]
prot2 <- aggregate(temp ~ Genes2, data = prot, FUN = mean)
prot <- prot2
colnames(prot) <- c('Genes','D492M/D492')
prot <- prot[-which(!is.finite(prot$`D492M/D492`)),]
prot$type <- 'prot'


# Now perform the microarray:
micro <- microarray[,c(1,6)]
colnames(micro)[2] <- 'temp'
micro <- micro[which(micro$PROBE_ID %in% hugo_genes$illumina_humanht_12_v3),]
micro$Genes2 <- NA
for (i in 1:nrow(micro)){
  if (micro$PROBE_ID[i] %in% hugo_genes$illumina_humanht_12_v3){
    micro$Genes2[i] <- unique(hugo_genes$hgnc_symbol[which(hugo_genes$illumina_humanht_12_v3 %in% micro$PROBE_ID[i])])
  }
}

micro <- data.frame(micro[,c(3,2)])
micro <- aggregate(temp ~ Genes2, data = micro, FUN = mean)
colnames(micro) <- c('Genes','D492M/D492')
micro$type <- 'micro'
list_names <- list(rnaseq$Genes,prot$Genes,micro$Genes)
list_of_genes <- Reduce(intersect,list_names)


# Find the common gene identifiers between datasets
rnaseq <- rnaseq[which(rnaseq$Genes %in% list_of_genes),]
prot <- prot[which(prot$Genes %in% list_of_genes),]
micro <- micro[which(micro$Genes %in% list_of_genes),]

# Combine all dataframes into a single dataframe:
df_fin <- rnaseq
colnames(df_fin)[2] <- 'rnaseq'
df_fin <- df_fin[,-3]
id_vec1 <- c()
for (i in 1:nrow(df_fin)){
  id_vec1 <- c(id_vec1,which(micro$Genes %in% df_fin$Genes[i]))
}
df_fin$micro <- NA
df_fin$micro <- micro$`D492M/D492`[id_vec1]

id_vec2 <- c()
for (i in 1:nrow(df_fin)){
  id_vec2 <- c(id_vec2,which(prot$Genes %in% df_fin$Genes[i]))
}
df_fin$prot <- NA
df_fin$prot <- prot$`D492M/D492`[id_vec2]
df_fin <- data.frame(df_fin)

df_fin$micro <- log(df_fin$micro,2)
df_fin$prot <- log(df_fin$prot,2)
df_fin$rnaseq <- log(df_fin$rnaseq,2)

colnames(df_fin)[c(2:4)] <- c('RNAseq','Microarray','Proteome')

#trace(chart.Correlation, edit=TRUE) # This was to change the correlation value sizes to be equal. 

#png(file = "Correlation_rnaseq_micro_prot_DEC2020.png",width = 6, height = 6, units = 'in', res = 300)
chart.Correlation(df_fin[,-1], histogram=TRUE, pch=20,method = 'spearman')
#dev.off()


# Now do the correlation analysis for metabolic identifiers only:

# Downloaded the metabolic gene list from REACTOME
metabolic <- read.csv("result (2).csv",header=T, sep=",")
subset <- metabolic[grep("^Metabolism$",metabolic$Pathway.name),]
metabolic_genes = subset$Mapped.entities
haha <- strsplit(as.character(metabolic_genes),";")
df <- data.frame(matrix(unlist(haha)))
metabolic_genes2 <- hugo_genes$hgnc_symbol[which(hugo_genes$uniprot_gn_id %in% df$matrix.unlist.haha..)]
df_fin_met <- df_fin[which(df_fin$Genes %in% metabolic_genes2),] # only about 300 genes here...


#png(file = "Correlation_rnaseq_micro_prot_MET_ONLY_DEC2020.png",width = 6, height = 6, units = 'in', res = 300)
chart.Correlation(df_fin_met[,-1], histogram=TRUE, pch=19,method = 'spearman')
#dev.off()