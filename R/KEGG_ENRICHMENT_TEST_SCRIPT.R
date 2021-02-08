# Proteomic data enrichment analysis:

# Upload the list of signifiacnatly differentially expressed proteins, and perform an
# enrichment analysis of the FDR-corrected p-values.

# Sigurdur Karvelsson

library(tidyr)
library(biomaRt)
library(readxl)
library(limma)
library(GO.db)
library(org.Hs.eg.db)
library(ggplot2)

prot <- read_excel('Analysis Results Qiong.xlsx')
prot <- prot[,c(2,5,10,11,13,14,16)]
prot <- data.frame(prot)
colnames(prot) <- c('Uniprot','EPI_1','MES_1','MES_2','EPI_2','EPI_3','MES_3')

# Take out proteins with more than 1 missing data in either condition:
prot <- prot[-which(apply(prot[,-1], 1, function(c)sum(c!=0)) <= 4),]
# Take out last value
prot <- prot[-nrow(prot),]

# Gather data
data_long <- gather(prot, condition, measurement,EPI_1:MES_3,factor_key=TRUE)
data_long$condition <- substring(data_long$condition,1,3)
data_long$condition <- as.factor(data_long$condition)
data_long$measurement[which(data_long$measurement == 0)] <- NA

glog <- function(x){
  res <- log((x + sqrt(x^2 + 1))/2,2)
  #res <- log((x + sqrt(x^2 + 1))/2,2)
  return(res)
}

# Log-transform data using generalized log transformation
data_long$measurement <- glog(data_long$measurement)
sum_dat <- data.frame(cbind(protein = prot$Uniprot,p_val =NA,adj_pval =NA))

# Perform series of t-tests:
for (i in 1:nrow(prot)){
  dat_new <- data_long[which(data_long$Uniprot %in% prot$Uniprot[i]),]
  temp_test <- t.test(measurement ~ condition,data = dat_new,na.action = na.omit)
  sum_dat$p_val[i] <- temp_test$p.value
}
sum_dat$p_val <- as.numeric(sum_dat$p_val)
sum_dat$adj_pval <- p.adjust(sum_dat$p_val,method = 'fdr',n = nrow(prot))

# Map the IDs:
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
hugo_genes <- getBM(attributes=c('entrezgene_id','uniprot_gn_id','illumina_humanht_12_v3'), 
                    filters = 'uniprot_gn_id', 
                    #values = colnames(rnaseq_D492_2), 
                    values = sum_dat$protein,
                    mart = ensembl)
sum_dat$gene_ids <- NA
for (i in 1:nrow(sum_dat)){
  if (sum_dat$protein[i] %in% hugo_genes$uniprot_gn_id){
    sum_dat$gene_ids[i] <- unique(hugo_genes$entrezgene_id[which(hugo_genes$uniprot_gn_id %in% sum_dat$protein[i])])
  } else {
    sum_dat$gene_ids[i] <- NA
  }
}

# PERFORM THE ENRICHMENT ANALYSIS!
Genes <- sum_dat$gene_ids[which(sum_dat$adj_pval < 0.05)]
g <- goana(Genes,species = 'Hs',universe = sum_dat$gene_ids)
topGO(g)

keg <- kegga(Genes, species="Hs",universe = sum_dat$gene_ids)
topKEGG(keg)

keg2 <- keg[which(keg$P.DE < 0.05),]
keg2 <- keg2[order(keg2$P.DE,decreasing = F),]
keg2$Pathway <- factor(keg2$Pathway, levels = rev(keg2$Pathway))

p <- ggplot(keg2, aes(x = Pathway, y = -log(P.DE,10),fill = DE/N))+
  geom_bar(stat = 'identity',width = 0.7) + scale_fill_gradient(low = "white", high = "grey30", 
                                                                na.value = NA,name = 'Enrichment factor')
p <- p + coord_flip()
p <- p + theme_bw() + ylab('-log(Enrichment p-value)') + ggtitle('KEGG - Enriched proteomic EMT pathways') + 
  theme(axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.title.x = element_text(size=12, face = "bold", colour = "grey10",
                                    margin = margin(t = 15, r = 15, b = 0, l = 0),vjust = 1),
        axis.title.y = element_blank(),
        axis.text.x = element_text( color="grey10", 
                                    size=10),
        axis.text.y = element_text(color="grey10", 
                                   size=12),
        legend.position = c(0.75,0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color='grey'),
        strip.text.x = element_text(size = 12))
p
ggsave('Enriched_KEGG_pathways_proteomic_data_JAN2021.png',width = 12, height = 5, units = "in")


