setwd("/home/oxygen/Ali_MB/")

library(Seurat) 
library(SeuratObject) 
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(tidydr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(DoubletFinder)

for (i in 1:4) {
  
  source("./R_function/ScRNAseq_preprocess_function_human.R")
  metadata <- readxl::read_xlsx("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/Clinical_metadata.xlsx") %>% 
    na.omit() %>%
    as.data.frame()
  
  
  gz_mat_path = paste0("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/",metadata$Sample_id[i])
  sample_id = metadata$Sample_id[i]
  doublet_rate = metadata$Doublet_rate[i]
  
  data <- readH5_QC_filter_doublet_remove(gz_mat_path = gz_mat_path,
                                          sample_id = sample_id,
                                          doublet_rate = doublet_rate)
  saveRDS(data, file = paste0("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/", metadata$Sample_id[i], ".rds"))
  
  rm(list = ls())
  
}

#integration
metadata <- readxl::read_xlsx("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/Clinical_metadata.xlsx") %>% 
  na.omit() %>%
  as.data.frame()

sc_list <- list()

for (i in 1:4) {
  sc_list[[i]] <- readRDS(paste0("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/",metadata$Sample_id[i],".rds"))
}

#####merge samples#####
features <- SelectIntegrationFeatures(object.list = sc_list)
anchors <- FindIntegrationAnchors(object.list = sc_list, 
                                  anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

# original unmodified data still resides in the 'RNA' assay
DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
# scale data
integrated <- ScaleData(integrated, verbose = FALSE)

# Add cell cycle score
integrated <- CellCycleScoring(integrated, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# Define variables in metadata to regress
vars_to_regress <- c("nFeature_RNA", "S.Score", "G2M.Score", "percent.mt", "percent.rb")

# Regress out the uninteresting sources of variation in the data
integrated <- ScaleData(object = integrated, vars.to.regress = vars_to_regress, verbose = T)

# clustering
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:15)
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:15)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:15)
integrated <- FindClusters(integrated, resolution = seq(0.1,2,0.1))
saveRDS(integrated, file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/integrated.rds")

#####annotation#####
library(Seurat) 
library(clustree)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(tidydr)

integrated <- readRDS("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/integrated.rds")

clustree_integrated <- clustree(integrated) + 
  guides(color = guide_legend(ncol = 3))
ggsave(clustree_integrated, width = 10, height = 10, filename = "./plots/QC/GSE141526_Human_Basal_cell_carcinoma_Clustree.pdf")

#optimal clustering resolution 0.3
integrated$seurat_clusters <- integrated$integrated_snn_res.0.2

theme_my <- function (xlength = 0.3, ylength = 0.3, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) {theme_classic() %+replace% theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = xlength, arrow = arrow), axis.line.y.left = element_line2(id = 2, ylength = ylength, arrow = arrow), axis.title = element_text(hjust = 0.1))}
cluster_color <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(length(unique(integrated$seurat_clusters)))

dimplot_cluster <- DimPlot(integrated, reduction = "tsne", label = T, repel = T, group.by = "seurat_clusters",
                           cols = cluster_color) + 
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))  +
  theme(legend.position = "none",
        title = element_blank())

dimplot_cluster <- DimPlot(integrated, reduction = "umap", label = T, repel = T, group.by = "seurat_clusters",
                           cols = cluster_color) + 
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))  +
  theme(legend.position = "none",
        title = element_blank())

dimplot_sample_id <- DimPlot(integrated, reduction = "tsne", group.by = "orig.ident",
                             cols = colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(length(unique(integrated$orig.ident)))) +
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(title = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, 'mm'))
p <- dimplot_cluster + dimplot_sample_id
ggsave(plot = p, filename = "./plots/QC/GSE141526_Human_Basal_cell_carcinoma_recluster_sample_id.pdf", width = 20, height = 9, units = "cm")

#identify clusters
features_set <- c(
  "KRT10", #Neoplastic cell
  "MLANA", "DCT",
  "CD79A","CD79B", # B cell
  "PTPRC", # immune
  "CD3E", "CD3D", "CD4", # CD4T cell
  "CD8A",# CD8T cell
  "VWF", "PLVAP",  # endothelial
  "COL1A1", "COL1A2", "LUM", "COL6A2", "FN1",#fibroblast
  "TPSAB1", #Mast cell
  "LYZ", "CD68", "TYROBP", # myeloid  lineage
  "RGS5", "ACTA2", "PDGFRB", # mural cell
  "GNLY", "NKG7", #NK cell
  "CSF3R" #Neutrophil
)

vlnplot_unannotated <- VlnPlot(integrated, assay = "RNA",
                               features = features_set, stack = T, flip = T,
                               group.by = "seurat_clusters",
                               cols = colorRampPalette(brewer.pal(n = 10, name = "Paired"))(length(features_set))) +
  NoLegend() +
  theme(strip.text.y = element_text(face = "italic", hjust = 0),
        axis.title.x = element_blank())
ggsave(plot = vlnplot_unannotated, filename = "./plots/QC/GSE141526_Human_Basal_cell_carcinoma_annotation_vlnplot_marker.pdf", width = 10, height = 12, units = "cm")

dotplot_unannotated <- DotPlot(integrated, features = rev(features_set), group.by = "seurat_clusters", dot.scale = 6, assay = "RNA") +
  coord_flip() +  
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.3) +
  viridis::scale_colour_viridis(option="magma", direction = 1) +
  guides(size=guide_legend("Perc. exp.", override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.y = element_text(face = "italic"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 30, size = 7))
ggsave(plot = dotplot_unannotated, filename = "./plots/QC/GSE141526_Human_Basal_cell_carcinoma_annotation_dotplot_marker.pdf", width = 14, height = 12, units = "cm")


p <- dimplot_cluster + vlnplot_unannotated
ggsave(plot = p, filename = "./plots/QC/GSE141526_Human_Basal_cell_carcinoma_annotation_dimplot_vlnplot_unanotated.pdf", width = 22, height = 12, units = "cm")

#annotate cluster
integrated$integrated_cell_type <- plyr::mapvalues(integrated$seurat_clusters,
                                                   from=0:10,
                                                   to=c("Keratinocyte",#0
                                                        "Keratinocyte",#1
                                                        "Fibroblast",#2
                                                        "Keratinocyte",#3
                                                        "Mural cell",#4
                                                        "Keratinocyte",#5
                                                        "Endothelial cell",#6
                                                        "Melanocyte",#7
                                                        "Keratinocyte",#8
                                                        "T cell",#9
                                                        "Keratinocyte"
                                                   ))


vlnplot <- VlnPlot(integrated, assay = "RNA",
                   features = features_set, stack = T, flip = T,
                   group.by = "integrated_cell_type",
                   cols = colorRampPalette(brewer.pal(n = 10, name = "Paired"))(length(features_set))) +
  NoLegend() +
  theme(strip.text.y = element_text(face = "italic", hjust = 0),
        axis.title.x = element_blank())
ggsave(plot = vlnplot, filename = "./plots/QC/GSE141526_Human_Basal_cell_carcinoma_vlnplot_marker.pdf", width = 10, height = 12, units = "cm")

dotplot <- DotPlot(integrated, features = rev(features_set), group.by = "integrated_cell_type", dot.scale = 6, assay = "RNA") +
  coord_flip() +  
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.3) +
  viridis::scale_colour_viridis(option="magma", direction = 1) +
  guides(size=guide_legend("Perc. exp.", override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.y = element_text(face = "italic"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(plot = dotplot, filename = "./plots/QC/GSE141526_Human_Basal_cell_carcinoma_dotplot_marker.pdf", width = 14, height = 12, units = "cm")


p <- dimplot_cluster + vlnplot + dotplot +
  dimplot_sample_id + vlnplot_unannotated + dotplot_unannotated +
  plot_layout(nrow = 2, ncol = 3)
ggsave(plot = p, filename = "./plots/QC/GSE141526_Human_Basal_cell_carcinoma_annotation_dimplot_vlnplot_anotated.pdf", width = 45, height = 28, units = "cm")

saveRDS(integrated, file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/integrated.rds")

#add SHH score to cells
#export expression matrix with 
library(clusterProfiler)
library(Seurat) 
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(tidydr)

cancer_hallmark <- read.gmt("/data/Shared_software/ref_genome/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt")
cancer_hallmark_table <- data.frame(pathway = unique(cancer_hallmark$term))

go_pathway <- read.gmt("/data/Shared_software/ref_genome/msigdb_v2023.1.Hs_GMTs/c5.all.v2023.1.Hs.symbols.gmt")
go_pathway_table <- data.frame(pathway = unique(go_pathway$term)) 

kegg_pathway <- read.gmt("/data/Shared_software/ref_genome/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
kegg_pathway_table <- data.frame(pathway = unique(kegg_pathway$term))


selected_pathway <- c(
  #Stress fiber
  "GOBP_STRESS_FIBER_ASSEMBLY", 
  "GOBP_POSITIVE_REGULATION_OF_STRESS_FIBER_ASSEMBLY", 
  
  #Mechanic stimulus response
  "GOBP_RESPONSE_TO_MECHANICAL_STIMULUS", 
  "GOBP_CELLULAR_RESPONSE_TO_MECHANICAL_STIMULUS", 
  
  #SHH pathway
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
  
  #proliferation
  "GOBP_STEM_CELL_PROLIFERATION",
  "GOBP_REGULATION_OF_STEM_CELL_PROLIFERATION"
  
)

#subset pathways
pathways <- do.call(rbind, list(
  cancer_hallmark %>% filter(term %in% selected_pathway),
  go_pathway %>% filter(term %in% selected_pathway),
  kegg_pathway %>% filter(term %in% selected_pathway),
  diapause,
  Stemness
))

# convert pathway to list
pathway_signature_list <- list()
for (pathway in unique(pathways$term)) {
  data <- pathways %>% 
    filter(term == pathway) 
  pathway_signature_list[[pathway]] <- data$gene
}


integrated[["joined"]] <- JoinLayers(integrated[["RNA"]])
DefaultAssay(integrated) <- "joined"

integrated  <- AddModuleScore(object = integrated , 
                              features = pathway_signature_list, 
                              name = names(pathway_signature_list),
                              assay = "joined")


# !!!! remove the tailing numbers of pathways added by seurat !!!!
# need to modify if order of pathway name changed !!!!
colnames(integrated@meta.data)[c(48:54)] <- gsub("\\d*$", "", colnames(integrated@meta.data)[c(48:54)])

saveRDS(integrated, file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/integrated.rds")
rm(list = ls())

#####dim plot, pathway plot#####
library(ComplexHeatmap) 
library(circlize)
library(Seurat) 
library(data.table)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(tidydr)

integrated <- readRDS("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/integrated.rds")

#plot dimplot feature plot
theme_my <- function(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) { theme_classic() %+replace% theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = xlength, arrow = arrow), axis.line.y.left = element_line2(id = 2, ylength = ylength, arrow = arrow), axis.title = element_text(hjust = 0.1))}
dimplot <- DimPlot(integrated, reduction = "tsne", group = "integrated_cell_type",
                   cols = colorRampPalette(brewer.pal(n = 11, name = "Paired"))(length(unique(integrated$integrated_cell_type))) )+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "") +
  theme_my()

mechano_response_plot <- FeaturePlot(integrated, features = c("GOBP_RESPONSE_TO_MECHANICAL_STIMULUS"), reduction = "tsne") +
  labs(title = "GOBP RESPONSE TO MECHANICAL STIMULUS") +
  scale_colour_gradient(low = "lightgray", high = "red2") +
  theme_my() + theme(title = element_text(size = 6))

design <- "AB"

p <- dimplot + mechano_response_plot +
  plot_layout(design = design)

ggsave(filename = "./plots/GSE141526_Human_Basal_cell_carcinoma_dimplot_mechano_response.pdf", plot = p, width = 7, height = 3)

##### pathway cor plot #####
# cor plot with * mark
# pval < 0.5
library(ComplexHeatmap) 
library(circlize)
library(Seurat) 
library(data.table)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(tidydr)

integrated <- readRDS("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/GSE141526_Human_Basal_cell_carcinoma_scRNA_Catherine_2020/integrated.rds")

#extract pathway score
pathway_expr <- integrated@meta.data
pathway_expr <- pathway_expr[,48:54]

#extract metadata
metadata <- integrated@meta.data

neoplastic_cell <- metadata %>%
  filter(integrated_cell_type %in% c("Keratinocyte"))
neoplastic_cell <- rownames(neoplastic_cell)

selected_pathway <- c(
  
  #Mechanic stimulus response
  "GOBP_RESPONSE_TO_MECHANICAL_STIMULUS", 
  
  #SHH pathway
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
  
  #Stress fiber
  "GOBP_STRESS_FIBER_ASSEMBLY", 
  "GOBP_POSITIVE_REGULATION_OF_STRESS_FIBER_ASSEMBLY", 
  
  #proliferation
  "GOBP_STEM_CELL_PROLIFERATION",
  "GOBP_REGULATION_OF_STEM_CELL_PROLIFERATION"
  
)

#arrange data
metadata_arrange <- metadata %>%
  filter(rownames(.) %in% neoplastic_cell) %>%
  arrange(desc(GOBP_RESPONSE_TO_MECHANICAL_STIMULUS)) 

pathway_expr <- pathway_expr[rownames(metadata_arrange), selected_pathway]

colnames(pathway_expr) <- c(
  #Mechanic stimulus response
  "GOBP Response to mechanical stimulus", 
  
  #SHH pathway
  "KEGG Hedgehog signaling pathway",
  
  #Stress fiber
  "GOBP Stress fiber assembly", 
  "GOBP Positive regulation of stress fiber assembly", 
  
  #proliferation
  "GOBP Stem cell proliferation",
  "GOBP Regulation of stem cell proliferation"
)

pathway_split <- c(
  rep("Mechano-response", 1),
  rep("Hedgehog signaling", 1),
  rep("Stress fiber", 2),
  rep("Stem cell proliferation", 2)
)

pathway_split <- factor(x = pathway_split, levels = c("Mechano-response", "Hedgehog signaling", "Stress fiber", "Stem cell proliferation"))

pathway_expr <- scale(pathway_expr)

pathway_expr_heatmap <- Heatmap(t(pathway_expr),
                                column_title=NULL,
                                row_split = pathway_split,
                                row_title_rot = 0,
                                show_column_dend = F,
                                show_column_names = F,
                                row_names_gp = gpar(fontsize = 6, col = "black"),
                                cluster_rows = F,
                                cluster_columns = F,
                                width = unit(10, "cm"), 
                                height = unit(5, "cm"),
                                name = "Pathway score",
                                border_gp = gpar(col = "black", lty = 1),
                                heatmap_legend_param = list(border = "black", direction = "horizontal",
                                                            legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), 
                                                            grid_width = unit(3, "mm")),
                                col = colorRamp2(c(-2, 0, 2), c("#fff500", "#040500", "#fe0000")),
                                raster_by_magick = FALSE,
                                use_raster = F
)

# cor heatmap
library(Hmisc)

# calculate correlation and p-values
res <- rcorr(as.matrix(pathway_expr), type = "pearson")

cor_mat <- res$r
p_mat <- res$P

# Replace NA p-values (diagonal) with 1
p_mat[is.na(p_mat)] <- 1

# cor with *
cor_heatmap <- Heatmap(cor_mat,
                       show_column_dend = F,
                       show_column_names = T,
                       column_split = pathway_split,
                       column_title = NULL,
                       cluster_rows = F,
                       cluster_columns = F,
                       width = unit(5, "cm"), 
                       height = unit(5, "cm"),
                       name = "PCC",
                       
                       # add asterisk 
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if(p_mat[i, j] < 0.05) {
                           grid.text("*", x, y, gp = gpar(fontsize = 15, col = "black"))
                         }
                       },
                       
                       border_gp = gpar(col = "black", lty = 1),
                       heatmap_legend_param = list(
                         border = "black", 
                         direction = "horizontal",
                         legend_height = unit(20, "mm"), 
                         labels_gp = gpar(fontsize = 8), 
                         grid_width = unit(3, "mm")
                       ),
                       col = colorRamp2(c(-1, 0, 1), c("#4c588a", "white", "#de2e29"))
)

heatmap_pathway_cor <- pathway_expr_heatmap + cor_heatmap
pdf(file = "./plots/GSE141526_Human_Basal_cell_carcinoma_heatmap_pathway_score_cor_with_pval_selected_pathway.pdf", width = 20, height = 20)
draw(heatmap_pathway_cor, heatmap_legend_side = "bottom")
dev.off()
