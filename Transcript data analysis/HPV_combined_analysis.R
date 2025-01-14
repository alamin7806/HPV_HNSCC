library(tidyverse)
library(VennDiagram)
library(gplots)
library(pheatmap)
library(RSkittleBrewer)
library(ggVennDiagram)
library(purrr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)

getwd()
setwd("C:/Users/hp/Desktop/HPV/GSE70462/RNAseq_data_analysis/GSE70463")
source("GSE70462.R")
setwd("C:/Users/hp/Desktop/HPV/GSE72536/RNASeq_data/GSE72536_rdata")
source("GSE72536.R")
setwd("C:/Users/hp/Desktop/HPV/GSE211322/RNA_data_analysi/RnA_seq322")
source("GSE211322.R")
setwd("C:/Users/hp/Desktop/HPV/GSE250305/RNASEQ_GSE250305")
source("GSE250305.R")

ensg462 = annotated_462df$ensgene
ensg536 = annotated_536df$ensgene
ensg322 = annotated_322df$ensgene
ensg305 = annotated_305df$ensgene

## intersect genomic data

common_HPV_ensg = Reduce(intersect, list(ensg305, ensg322, ensg462, ensg536))

ensg462_commonGene = annotated_462df[annotated_462df$ensgene
      %in% common_HPV_ensg, c("ensgene", "log2FoldChange")]

ensg536_commonGene = annotated_536df[annotated_536df$ensgene
                        %in% common_HPV_ensg,c("ensgene", "log2FoldChange")]

ensg322_commonGene = annotated_322df[annotated_322df$ensgene
                       %in% common_HPV_ensg,c("ensgene", "log2FoldChange")]
ensg305_commonGene = annotated_305df[annotated_305df$ensgene
                      %in% common_HPV_ensg,c("ensgene", "log2FoldChange")]

list_df = list(ensg305_commonGene, ensg322_commonGene, ensg462_commonGene,
               ensg536_commonGene)




common_folChange = purrr::reduce(list_df, left_join , by = "ensgene")
colnames(common_folChange)[c(2,3,4,5)] = c("GSE250305_log2foldchange", 
                                           "GSE211322_log2foldchange",
                                           "GSE70462_log2foldchange",
                                           "GSE72536_log2foldchange")

head(common_folChange)
view(list_df)
correlate_gene = cor(common_folChange[,c("GSE250305_log2foldchange", 
                                         "GSE211322_log2foldchange",
                                         "GSE70462_log2foldchange",
                                         "GSE72536_log2foldchange")], 
                     method = "pearson")


correlate_gene
view(correlate_gene)

anno305_diff = anno305_df3$ensgene
anno322_diff = anno322_df3$ensgene
anno462_diff = anno462_df3$ensgene
anno536_diff = anno536_df3$ensgene

## Common significant genes

common_gene_diff = Reduce(intersect, list(anno305_diff, anno322_diff,
                      anno462_diff, anno536_diff))

list_DE_gene = list(GSE250305= anno305_diff, 
                    GSE211322= anno322_diff, 
                    GSE70462= anno462_diff, 
                    GSE72536= anno536_diff)

v = venn(list_DE_gene)

v2 = ggVennDiagram(list_DE_gene) +
ggtitle("Venn diagram of overlapping DEGs frominter section of four independent GEO datasets")
v2

union_DE_genes = Reduce(union, list_DE_gene)

union_foldcchange_hm = filter(common_folChange, ensgene %in% union_DE_genes)
union_foldcchange_hm_mat = as.matrix(union_foldcchange_hm[,2:5])

pheatmap(union_foldcchange_hm_mat)

colramp = colorRampPalette(c(3,"white", 2))(9)

my_breaks = c(seq(-2, -0.01, length.out=50),
              0,
              seq(0.01, 2, length.out=50))


pheatmap(union_foldcchange_hm_mat, breaks = my_breaks, color = colramp,
         main = "Overlapping Genes from different GEO dataset")


GO_results <- enrichGO(gene = common_gene_diff,
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "ENSEMBL", 
                       ont = "BP")

entgene_DE = getBM(attributes = c("entrezgene_id"),
                    filters =c("ensembl_gene_id"),
                    values = common_gene_diff,
                    mart = ensemble111)

entgene_DE = as.character(entgene_DE$entrezgene_id)

entUni_DE = getBM(attributes = c("entrezgene_id"),
                   filters =c("ensembl_gene_id"),
                   values = common_HPV_ensg,
                   mart = ensemble111)
entUni_DE = as.character(entUni_DE$entrezgene_id)

ego_DE = enrichGO(gene = entgene_DE,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   universe = entUni_DE,
                   readable = TRUE)

ego_DE
view(summary(ego_DE))
barplot(ego_DE, title = "Overlapping Genes")
dotplot(ego_DE, title = "Overlapping Genes")




ekegg_DE = enrichKEGG(gene = entgene_DE,
                      universe = entUni_DE)

view(ekegg_DE)

## individual gene expression

plotCounts(dds322, gene = "ENSG00000178222", 
           intgroup = "condition", main = "GSE211322")

plotCounts(dds305, gene = "ENSG00000178222", 
           intgroup = "condition", main = "GSE250305")

plotCounts(dds462, gene = "ENSG00000178222", 
           intgroup = "condition", main = "GSE70462")

plotCounts(dds536, gene = "ENSG00000178222", 
           intgroup = "condition", main = "GSE72536")





















































