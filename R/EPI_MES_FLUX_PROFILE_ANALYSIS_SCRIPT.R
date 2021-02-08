# Script to perform two things: 1) Analyse the correlation in flux profiles between the different data type-constrained
# GSMMs. Subsequently, plot them up in a heatmap. 2) Analyse the compartment-related correlation in flux differences 
# between D492 and D492M GSMMs. 

# Sigurdur Karvelsson

library(readxl)
library(siggitRausti)
library(data.table)
library(tidyr)
library(ggpubr)
library(Hmisc)
library(ComplexHeatmap)



# 1) Plot the median flux profiles (from random sampling) of all GSMMs in a heatmap:
####################################################################################
labels <- read.table('median_rxns_all_GSMMs_compartment_DEC2020.txt',header=T,sep='\t')
colnames(labels) = labels[1,]
labels = labels[-1,]
for (i in 2:ncol(labels)){
  labels[,i] <- as.numeric(labels[,i])
}

cytosol = labels[which(labels$Location %in% '[c]'),]
labels <- labels[-which(labels$Location == 'transport'),]
colnames(labels)[2:ncol(labels)] <- c('Proteomic-EPI','Proteomic-MES','Microarray-EPI','Microarray-MES','Media-EPI','Media-MES',
                                      'RNAseq-EPI','RNAseq-MES')

# Create a heatmap of all reactions and annotate with subsystems:
col = list(Phenotype = c("EPI" = "dodgerblue4", "MES" = "firebrick4"))

# Make a new labels dataframe with scaled values:
labels2 <- data.frame(t(labels[,-1]))
labels2 <- data.frame(scale(labels2))

text_list = list(
  text1 = 'Cytosol',
  text2 = 'Mitochondria',
  text3 = 'ER',
  text4 = 'Lysosome',
  text5 = 'Golgi',
  text6 = 'Peroxisome',
  text7 = 'Nucleus',
  text8 = 'Extracellular'
)

ha = rowAnnotation(foo = anno_empty(border = F, 
                                    width = max_text_width(unlist(text_list)) + unit(4, "mm")))


ha2 <- HeatmapAnnotation(
  Phenotype = c('EPI','MES','EPI','MES','EPI','MES','EPI','MES'),
  annotation_label = c('Phenotype'),
  col = col,
  gp = gpar(col = "black"),
  border = T
)

heatmap2 <- Heatmap(as.matrix(t(labels2)), name = "Reaction activity",
                    row_split = factor(labels$Location,levels = c('[c]','[m]','[r]','[l]','[g]','[x]','[n]','[e]')),
                    row_title = ' ',
                    cluster_row_slices = F,
                    column_names_rot = 45,
                    row_names_gp = gpar(fontsize = 8),
                    right_annotation = ha,top_annotation = ha2,
                    show_row_names = F,
                    clustering_method_columns = 'ward.D2',
                    border=T)
heatmap2
#png(file = "compartment_heatmap_DEC2020.png",width = 6, height = 8, units = 'in', res = 300)
heatmap2
draw(heatmap2, auto_adjust = FALSE,merge_legend = T)
for(i in 1:8) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
  })
}
#dev.off()


# 2) Compare the correlation between median EPI and MES flux profiles between different
# data types. Check if there is a difference in the correlation of MES/EPI flux between 
# different cellular compartments. 
#######################################################################################
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# Now generate a heatmap of correlated reactions based on their location within the cells. 
labels_relative <- labels
labels_relative$prot <- labels_relative$`Proteomic-MES`/labels_relative$`Proteomic-EPI`
labels_relative$micro <- labels_relative$`Microarray-MES`/labels_relative$`Microarray-EPI`
labels_relative$media <- labels_relative$`Media-MES`/labels_relative$`Media-EPI`
labels_relative$rnaseq <- labels_relative$`RNAseq-MES`/labels_relative$`RNAseq-EPI`

compartments <- c('[c]','[m]','[r]','[l]','[g]','[x]','[n]','[e]')
compartment_cor_mat <- data.frame(cbind(Compartments = compartments,Average_spearman_cor = NA,pval_spearman_cor = NA))

# Check specifically the correlation between proteomic and RNAseq GSMMs:
for (i in 1:length(compartments)){
  comp_temp <- labels_relative[which(labels_relative$Location == compartments[i]),c(10,13)]
  res2 <- rcorr(as.matrix(comp_temp),type = 'spearman')
  res3 <- flattenCorrMatrix(res2$r, res2$P)
  compartment_cor_mat$Average_spearman_cor[i] <- res3$cor
  #compartment_cor_mat$std_spearman_cor[i] <- sd(res3$cor)
  compartment_cor_mat$pval_spearman_cor[i] <- res3$p
}
compartment_cor_mat$pval_spearman_cor <- as.numeric(compartment_cor_mat$pval_spearman_cor)
compartment_cor_mat$Average_spearman_cor <- as.numeric(compartment_cor_mat$Average_spearman_cor)
compartment_cor_mat$adj_pval_spearman_cor <- p.adjust(compartment_cor_mat$pval_spearman_cor,method = 'bonferroni')

compartment_cor_mat$pval_spearman_cor <- as.character(signif(compartment_cor_mat$pval_spearman_cor,2))
compartment_cor_mat$adj_pval_spearman_cor <- as.character(signif(compartment_cor_mat$adj_pval_spearman_cor,2))
compartment_cor_mat$Average_spearman_cor <- as.character(signif(compartment_cor_mat$Average_spearman_cor,2))
compartment_cor_mat
#write.table(compartment_cor_mat[,-1],'compartment_correlations_rnaseq_proteomic.txt',row.names = F,sep = '\t')