setwd("/home/oxygen/Ali_MB/")

library(tidyverse)
library(Rtsne)

#read expr
expr <- data.table::fread("/data/Sequencing_data_inhouse/Medulloblastoma/GSE124814_Human_MB_BulkRNA_Weishaupt_2019/GSE124814/GSE124814_HW_expr_matrix.tsv.gz") %>%
  as.data.frame() %>%
  column_to_rownames("Gene_Symbol")
metadata <- readxl::read_xlsx("/data/Sequencing_data_inhouse/Medulloblastoma/GSE124814_Human_MB_BulkRNA_Weishaupt_2019/GSE124814/GSE124814_sample_descriptions.xlsx", skip = 1) %>%
  as.data.frame() 

#select top variation gene
gene_sd <- data.frame(
  gene = rownames(expr),
  sd <- apply(expr, 1, sd)
) %>% slice_max(order_by = sd, n = 1500)
  
#pca
data <- t(expr[gene_sd$gene,]) # subset expr
pca <- prcomp(data) # row represents samples
pca_summary <- summary(pca)
pca_summary <- pca_summary$importance
plot(pca_summary[2,1:20])

#tsne
set.seed(2333)
tsne_xy <- Rtsne(
  data,
  dims = 2,
  initial_dims = 6,
  pca = T,
  max_iter = 1000,
  theta = 0.5,
  perplexity = 100,
  num_threads = 64,
  verbose = T
)  

tsne_xy <- tsne_xy$Y
colnames(tsne_xy) <- c("tSNE_1", "tSNE_2")

#add tsne to metadata
metadata <- cbind(metadata, tsne_xy)

#save processed data
metadata$subgroup <- ifelse(metadata$`source name` == "Normal", "Normal", metadata$`characteristics: subgroup relabeled`)
save(metadata, file = "/data/Sequencing_data_inhouse/Medulloblastoma/GSE124814_Human_MB_BulkRNA_Weishaupt_2019/GSE124814/metadata.Rdata")
save(expr, file = "/data/Sequencing_data_inhouse/Medulloblastoma/GSE124814_Human_MB_BulkRNA_Weishaupt_2019/GSE124814/expr_quantile_normalized.Rdata")

##### expression heatmap #####
library(ComplexHeatmap)
library(circlize)

# mechanosenser
Mechanosensor <- c(
  "PIEZO1",
  "PIEZO2",
  "TMEM63A",
  "TMEM63B",
  "SCNN1A",
  "TRPC1",
  "TRPC6",
  "TRPM7",
  "TRPV4",
  "PKD2L1",
  "KCNK2",
  "KCNK10",
  "TMC1",
  "TMC2",
  "GPR68"
)

Mechanosensor_mean_expr <- aggregate(t(expr[Mechanosensor,]), by = list(metadata$subgroup), FUN = mean) %>%
  column_to_rownames("Group.1")

max <- max(Mechanosensor_mean_expr)
min <- min(Mechanosensor_mean_expr)
middle <- (max + min) / 2

Mechanosensor_mean_expr <- Mechanosensor_mean_expr[c("SHH", "WNT", "G3", "G4", "Normal"), c(
  "GPR68",
  "TMEM63A",
  "PIEZO2",
  "TRPM7",
  "TRPC6",
  "TRPC1",
  "PKD2L1",
  "TMC1",
  "SCNN1A",
  "TRPV4",
  "KCNK2",
  "TMC2",
  "TMEM63B",
  "PIEZO1"
  )]
export::table2excel(Mechanosensor_mean_expr, "./processed_data/MB_bulk_Mechanosensor_mean_expr_by_subtype.xlsx")

Sample_count <- table(metadata$subgroup) %>% as.data.frame() %>% column_to_rownames("Var1")
Sample_count <- Sample_count[c("SHH", "WNT", "G3", "G4", "Normal"),]
right_annotation <- rowAnnotation(`Sample count` = anno_barplot(Sample_count, add_numbers = T, numbers_rot = 0, axis = F, gp = gpar(fill = c("#ea4335", "#673ab7", "#4285f4", "#fbbc05", "#34a853"))), 
                                  annotation_name_rot = 0,
                                  annotation_name_gp = gpar(fontsize = 7),
                                  width = unit(1.8, "cm"))

heatmap_final_version <- Heatmap(Mechanosensor_mean_expr,
                   column_split = 1:ncol(Mechanosensor_mean_expr),
                   column_title = colnames(Mechanosensor_mean_expr),
                   column_title_rot = 30,
                   column_title_gp = gpar(fontsize = 9, fontface = c(rep("bold.italic", 1), rep("italic", 13)), col = c(rep("red", 1), rep("black", 13))),
                   row_split = 1:nrow(Mechanosensor_mean_expr),
                   row_title = rownames(Mechanosensor_mean_expr),
                   row_title_rot = 0,
                   row_title_gp = gpar(fontsize = 11),
                   show_column_dend = F,
                   show_column_names = F,
                   show_row_names = F,
                   cluster_rows = F,
                   cluster_columns = F,
                   width = unit(14 * 0.75, "cm"), 
                   height = unit(5 * 0.5, "cm"),
                   right_annotation = right_annotation,
                   name = "Mean\nexpression",
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(ifelse(Mechanosensor_mean_expr[i, j] == max(Mechanosensor_mean_expr[,j]), "*", ""), x, y * 0.35, gp = gpar(fontsize = 20))
                   },
                   rect_gp = gpar(col = "black", lwd = 1),
                   col = colorRamp2(c(min, middle, max), c("white", "#efefef", "#af3026"))
)

##### correlation heatmap Mechanosenser/SHH gene in all tumor sample#####
library(clusterProfiler)

kegg_SHH <- read.gmt("/data/Shared_software/ref_genome/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt") %>%
  filter(term == "KEGG_HEDGEHOG_SIGNALING_PATHWAY") 

#subset gene
kegg_SHH <- kegg_SHH$gene[kegg_SHH$gene %in% rownames(expr)]

#make empty dataframe 
cor_SHH_mechanosensor <- data.frame(matrix(nrow = length(kegg_SHH), ncol = length(Mechanosensor)))
rownames(cor_SHH_mechanosensor) <- kegg_SHH
colnames(cor_SHH_mechanosensor) <- Mechanosensor

pval_SHH_mechanosensor <- data.frame(matrix(nrow = length(kegg_SHH), ncol = length(Mechanosensor)))
rownames(pval_SHH_mechanosensor) <- kegg_SHH
colnames(pval_SHH_mechanosensor) <- Mechanosensor

# select all tumor sample
metadata_tumor <- metadata %>% filter(`source name` == "Medulloblastoma")

# calculate cor and pval
for (row in kegg_SHH) {
  for (column in Mechanosensor) {
    cor_SHH_mechanosensor[row, column] <- cor.test(unlist(expr[row,metadata_tumor$`Sample name`]), unlist(expr[column,metadata_tumor$`Sample name`]), method = "pearson")[["estimate"]][["cor"]]
    pval_SHH_mechanosensor[row, column] <- cor.test(unlist(expr[row,metadata_tumor$`Sample name`]), unlist(expr[column,metadata_tumor$`Sample name`]), method = "pearson")[["p.value"]]
  }
}

export::table2excel(cor_SHH_mechanosensor, "./processed_data/MB_bulk_Mechanosensor_SHH_pearson_cor_coefficient_in_all_tumor_sample.xlsx")
export::table2excel(pval_SHH_mechanosensor, "./processed_data/MB_bulk_Mechanosensor_SHH_pearson_cor_pval_in_all_tumor_sample.xlsx")

cor_SHH_mechanosensor <- cor_SHH_mechanosensor %>%
  arrange(desc(GPR68))
pval_SHH_mechanosensor <- pval_SHH_mechanosensor[rownames(cor_SHH_mechanosensor),]

# make heatmap
max <- max(cor_SHH_mechanosensor)
min <- min(cor_SHH_mechanosensor)
middle <- (max + min) / 2

left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "#a2cf6e"),
                                                  labels = c("Correlation between mechanosenser and SHH pathway genes in all tumor"), 
                                                  labels_gp = gpar(col = "black", fontsize = 10)),
                                 width = unit(0.6, units = "cm")
)

heatmap2_final_version <- Heatmap(cor_SHH_mechanosensor[c("GLI1", "HHIP", "PTCH1", "GLI2", "PTCH2", "GLI3", "SMO", "GAS1"),],
                    column_split = 1:ncol(cor_SHH_mechanosensor),
                    column_title = rep("", ncol(cor_SHH_mechanosensor)),
                    #column_title_rot = 30,
                    #column_title_gp = gpar(fontsize = 9, fontface = c(rep("italic", 6), rep("bold.italic", 4), rep("italic", 4)), col = c(rep("black", 6), rep("red", 4), rep("black", 4))),
                    #row_split = 1:nrow(cor_SHH_mechanosensor),
                    #row_title = rownames(cor_SHH_mechanosensor),
                    #row_title_rot = 0,
                    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                    show_column_dend = F,
                    show_column_names = F,
                    show_row_names = T,
                    row_names_side = "left",
                    #left_annotation = left_annotation,
                    cluster_rows = F,
                    cluster_columns = F,
                    width = unit(14 * 0.75, "cm"), 
                    height = unit(3, "cm"),
                    #right_annotation = right_annotation,
                    name = "Pearson\ncorrelation\ncoefficient",
                    #cell_fun = function(j, i, x, y, width, height, fill) {
                    #  grid.text(ifelse( (abs(cor_SHH_mechanosensor[i, j]) > 0.5) & (pval_SHH_mechanosensor[i,j] < 0.05), sprintf("%.2f", cor_SHH_mechanosensor[i, j]), ""), x, y, gp = gpar(fontsize = 6))
                    #},
                    rect_gp = gpar(col = "black", lwd = 0.5),
                    col = colorRamp2(c(min, middle, max), c("#4c588a", "white", "#de2e29"))
)
heatmap2_final_version

##### correlation heatmap Mechanosenser/SHH pathway score in all tumor sample#####
# calculate SHH signature score
library(GSVA)
gsva_res <- gsva(expr = as.matrix(expr), gset.idx.list = list(SHH = kegg_SHH), method='ssgsea',kcdf='Gaussian')

gsva_res <- ssgseaParam(exprData = as.matrix(expr), geneSets = list(SHH = kegg_SHH))
gsva_res <- gsva(gsva_res)

#make empty dataframe 
cor_SHH_score_mechanosensor <- data.frame(matrix(nrow = 1, ncol = length(Mechanosensor)))
rownames(cor_SHH_score_mechanosensor) <- "SHH pathway score(all subgroups)"
colnames(cor_SHH_score_mechanosensor) <- Mechanosensor

pval_SHH_score_mechanosensor <- data.frame(matrix(nrow = 1, ncol = length(Mechanosensor)))
rownames(pval_SHH_score_mechanosensor) <- "SHH pathway score(all subgroups)"
colnames(pval_SHH_score_mechanosensor) <- Mechanosensor

# calculate cor and pval
for (row in 1) {
  for (column in Mechanosensor) {
    cor_SHH_score_mechanosensor[1, column] <- cor.test(unlist(gsva_res[1,metadata_tumor$`Sample name`]), unlist(expr[column,metadata_tumor$`Sample name`]), method = "pearson")[["estimate"]][["cor"]]
    pval_SHH_score_mechanosensor[1, column] <- cor.test(unlist(gsva_res[1,metadata_tumor$`Sample name`]), unlist(expr[column,metadata_tumor$`Sample name`]), method = "pearson")[["p.value"]]
  }
}

heatmap3_final_version <- Heatmap(cor_SHH_score_mechanosensor,
                    column_split = 1:ncol(cor_SHH_score_mechanosensor),
                    column_title = NULL,
                    row_title = "Correlation with\nSHH pathway",
                    row_title_rot = 0,
                    row_title_gp = gpar(fontsize = 8),
                    show_column_dend = F,
                    show_column_names = F,
                    show_row_names = F,
                    cluster_rows = F,
                    cluster_columns = F,
                    width = unit(14 * 0.75, "cm"), 
                    height = unit(0.5, "cm"),
                    show_heatmap_legend = F,
                    rect_gp = gpar(col = "black", lwd = 1),
                    col = colorRamp2(c(min, middle, max), c("#4c588a", "white", "#de2e29"))
)

##### export heatmap with all SHH pathway gene#####
ht_list = heatmap_final_version %v% heatmap3_final_version %v% heatmap2_final_version

pdf(file = "./plots/Bulk_MB_mechanosensor_SHH_heatmap_selected_gene_clean_version.pdf", width = 10, height = 8)
draw(ht_list, ht_gap = unit(c(0, 0.1, 0.1), "cm"))
dev.off()
