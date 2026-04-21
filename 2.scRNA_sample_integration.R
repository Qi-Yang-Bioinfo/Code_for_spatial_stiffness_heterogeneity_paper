setwd("/home/oxygen/Ali_MB/")
library(Seurat) 
library(SeuratObject) 
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(tidydr)
library(clustree)
library(BiocParallel)
library(data.table)

#####add metadata#####
#load samples
load(file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Experiment1/WT_core/preprocessed_doublet_filtered_WT_core_Experiment1.Rdata")
load(file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Experiment1/WT_peripheral_adjacent_to_skull/preprocessed_doublet_filtered_WT_peripheral_adjacent_to_skull_Experiment1.Rdata")

load(file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Experiment2/Gpr68KO_core/preprocessed_doublet_filtered_Gpr68KO_core_Experiment2.Rdata")
load(file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Experiment2/Gpr68KO_peripheral/preprocessed_doublet_filtered_Gpr68KO_peripheral_Experiment2.Rdata")

# add metadata
WT_core_Experiment1$batch <- "Experiment1"
WT_peripheral_adjacent_to_skull_Experiment1$batch <- "Experiment1"
Gpr68KO_core_Experiment2$batch <- "Experiment2"
Gpr68KO_peripheral_Experiment2$batch <- "Experiment2"

WT_core_Experiment1$genotype <- "WT"
WT_peripheral_adjacent_to_skull_Experiment1$genotype <- "WT"
Gpr68KO_core_Experiment2$genotype <- "Gpr68KO"
Gpr68KO_peripheral_Experiment2$genotype <- "Gpr68KO"
  
WT_core_Experiment1$sample_region <- "Core"
WT_peripheral_adjacent_to_skull_Experiment1$sample_region <- "Peripheral"
Gpr68KO_core_Experiment2$sample_region <- "Core"
Gpr68KO_peripheral_Experiment2$sample_region <- "Peripheral"

#####merge samples#####
features <- SelectIntegrationFeatures(object.list = c(
  WT_core_Experiment1,
  WT_peripheral_adjacent_to_skull_Experiment1,
  Gpr68KO_core_Experiment2,
  Gpr68KO_peripheral_Experiment2
))
anchors <- FindIntegrationAnchors(object.list = c(
  WT_core_Experiment1,
  WT_peripheral_adjacent_to_skull_Experiment1,
  Gpr68KO_core_Experiment2,
  Gpr68KO_peripheral_Experiment2
), anchor.features = features)
MB_integrated <- IntegrateData(anchorset = anchors)

# Run the standard workflow for visualization and clustering
# scale data
MB_integrated <- ScaleData(MB_integrated, verbose = FALSE)

# Add cell cycle score
hs_mm_dm_ortholog <- data.table::fread("./Human_orthologs/HS_symbol_MM_ortholog_DM_ortholog_ensembl_v110.txt") %>% as.data.frame()
mm.s.genes <- hs_mm_dm_ortholog %>% filter(`Gene name` %in% cc.genes$s.genes) %>% filter(`Mouse gene name` != "")
mm.s.genes <- mm.s.genes$`Mouse gene name`
mm.g2m.genes <- hs_mm_dm_ortholog %>% filter(`Gene name` %in% cc.genes$g2m.genes) %>% filter(`Mouse gene name` != "")
mm.g2m.genes <- mm.g2m.genes$`Mouse gene name`

DefaultAssay(MB_integrated) <- "RNA"
MB_integrated <- CellCycleScoring(MB_integrated, s.features = mm.s.genes, g2m.features = mm.g2m.genes)

# original unmodified data still resides in the 'RNA' assay
DefaultAssay(MB_integrated) <- "integrated"

# Define variables in metadata to regress
vars_to_regress <- c("nFeature_RNA", "S.Score", "G2M.Score", "percent.rb", "percent.mt")

# Regress out the uninteresting sources of variation in the data
MB_integrated <- ScaleData(object = MB_integrated, vars.to.regress = vars_to_regress, verbose = T)

# clustering
MB_integrated <- RunPCA(MB_integrated, npcs = 30, verbose = FALSE)
ElbowPlot(MB_integrated, ndims = 30, reduction = "pca")
ggsave(filename = "./plots/QC/Elbowplot.pdf", width = 5, height = 4)
MB_integrated <- RunUMAP(MB_integrated, reduction = "pca", dims = 1:14)
MB_integrated <- RunTSNE(MB_integrated, reduction = "pca", dims = 1:14)
MB_integrated <- FindNeighbors(MB_integrated, reduction = "pca", dims = 1:14)
MB_integrated <- FindClusters(MB_integrated, resolution = seq(0.1,2,0.1))

# export ncount ngene table per sample after integration
export::table2excel(do.call("cbind", tapply(MB_integrated$nCount_RNA, MB_integrated$orig.ident,quantile,probs=seq(0,1,0.05))), file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Sample_nCount_quantile.xlsx")
export::table2excel(do.call("cbind", tapply(MB_integrated$nFeature_RNA, MB_integrated$orig.ident,quantile,probs=seq(0,1,0.05))), file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Sample_ngene_quantile.xlsx")

#####Plot after integration, calculate markers for clusters#####
library(Seurat) 
library(SeuratObject) 
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(tidydr)
library(clustree)
library(BiocParallel)
library(data.table)
library(ggtext)

load("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/MB_integrated.Rdata")

# cluster tree plot
clustree_MB_integrated <- clustree(MB_integrated) + 
  scale_color_manual(values = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(20)) + 
  guides(color = guide_legend(ncol = 3))
ggsave(width = 20, height = 20, filename = "./plots/MB_integrated_Clustree.pdf")
# optimal resolution 0.8

# set seurat cluster
MB_integrated$seurat_clusters <- MB_integrated$integrated_snn_res.0.8

#plot previous annotation
theme_my <- function (xlength = 0.3, ylength = 0.3, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) 
{
  theme_classic() %+replace% theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = xlength, arrow = arrow), axis.line.y.left = element_line2(id = 2,
                                                                                                                                                         ylength = ylength, arrow = arrow), axis.title = element_text(hjust = 0.1))
}

orig.ident_dimplot <- DimPlot(MB_integrated, reduction = "tsne", group.by = "orig.ident", raster.dpi = c(1000, 1000), pt.size = 0.2,
                             cols = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(MB_integrated$orig.ident))) ) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 2.5), ncol=1)) + #smaller legend size
  theme( legend.text = element_text(size = 5))
orig.ident_dimplot

location_dimplot <- DimPlot(MB_integrated, reduction = "tsne", group.by = "sample_region", raster.dpi = c(1000, 1000), pt.size = 0.2,
                            cols = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(MB_integrated$sample_region))) ) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 2.5), ncol=1)) + #smaller legend size
  theme( legend.text = element_text(size = 5))
location_dimplot

genotype_dimplot <- DimPlot(MB_integrated, reduction = "tsne", group.by = "genotype", raster.dpi = c(1000, 1000), pt.size = 0.2,
                            cols = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(MB_integrated$genotype))) ) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 2.5), ncol=1)) + #smaller legend size
  theme( legend.text = element_text(size = 5))
genotype_dimplot

cluster_dimplot <- DimPlot(MB_integrated, reduction = "tsne", group.by = "integrated_snn_res.0.8", raster.dpi = c(1000, 1000), pt.size = 0.2,
                           cols = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(MB_integrated$integrated_snn_res.0.8))) ) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 2.5), ncol=1)) + #smaller legend size
  theme( legend.text = element_text(size = 5))
cluster_dimplot

p = orig.ident_dimplot + cluster_dimplot +
  location_dimplot + genotype_dimplot +
  plot_layout(nrow = 2, ncol = 2)

ggsave(filename = "./plots/MB_integrated_orig.ident_cluster_location_genotype.png", width = 13, height = 10, plot = p)

#calculate marker gene for clusters
markers <- FindAllMarkers(MB_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.7")
save(markers, file = "./processed_data/MB_integrated_snn_res.0.8_markers.Rdata")
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
export::table2excel(top_markers, file = "./processed_data/MB_integrated_snn_res.0.8_top20_markers.xlsx")

#####Add SHH A B C signature score#####
# calculate SHH A B C subtype
load("./external_data/GSE119926_Human_MB_scRNA_Hovestadt_2019_transcriptional_program.Rdata")
signature_list <- signature_list[c("SHH-A", "SHH-B", "SHH-C")]

# translate HS gene to MM gene
hs_mm_dm_ortholog <- data.table::fread("/data/Shared_software/ref_genome/drosophila/Human_orthologs/HS_symbol_MM_ortholog_DM_ortholog_ensembl_v110.txt") %>% as.data.frame()
signature_list_mm <- list()
for (signature in names(signature_list)) {
  gene_list <- signature_list[[signature]] 
  gene_list_mm <- hs_mm_dm_ortholog$`Mouse gene name`[match(gene_list, hs_mm_dm_ortholog$`Gene name`)]
  gene_list_mm[which(gene_list_mm == "")] <- NA
  signature_list_mm[[signature]] <- gene_list_mm[!is.na(gene_list_mm)]
}

MB_integrated <- AddModuleScore(object = MB_integrated, features = signature_list_mm, 
                                name = names(signature_list_mm), assay = "RNA")

cell_type <- MB_integrated@meta.data[, c("Sample", "Methylation subgroup", "WNT-A1", "WNT-B2", "WNT-C3", "WNT-D4", "SHH-A5", "SHH-B6", "SHH-C7", "Group 3/4-A8", "Group 3/4-B9", "Group 3/4-C10")]
cell_type$subtype <- NA

for (row in 1:nrow(cell_type)) {

  if (cell_type[row, "Methylation subgroup"] == "SHH") {
    max <- max(cell_type[row, c("SHH-A1", "SHH-B2", "SHH-C3")])
    cell_type$subtype[row] <- c("SHH-A", "SHH-B", "SHH-C")[which(cell_type[row, c("SHH-A5", "SHH-B6", "SHH-C7")] == max)]
  }
}

##### add signature score#####
library(clusterProfiler)

#read pathway
cancer_hallmark <- read.gmt("./msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt")
cancer_hallmark_table <- data.frame(pathway = unique(cancer_hallmark$term))

go_pathway <- read.gmt("./msigdb_v2023.1.Hs_GMTs/c5.all.v2023.1.Hs.symbols.gmt")
go_pathway_table <- data.frame(pathway = unique(go_pathway$term)) 

kegg_pathway <- read.gmt("./msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
kegg_pathway_table <- data.frame(pathway = unique(kegg_pathway$term))
  
Stemness <- readxl::read_xlsx("./external_data/Stemness_Alex_PNAS_2019.xlsx", sheet = 1) %>%
  as.data.frame() %>%
  mutate(term = "Stemness") %>%
  rename(gene = `Curated  (withouth immune and proliferative genes, annotated in TCGA)`) %>%
  select(term, gene) %>%
  na.omit()

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
  Stemness
))

# translate HS gene to MM gene
pathways$mm_gene <- hs_mm_dm_ortholog$`Mouse gene name`[match(pathways$gene, hs_mm_dm_ortholog$`Gene name`)]
pathways$na <- ifelse(pathways$mm_gene == "", 1, 0)
pathways$na <- ifelse(is.na(pathways$mm_gene), 1, pathways$na)
pathways$count <- 1

pathways_statistic <- aggregate(pathways$na, by = list(pathways$term), FUN = sum)
pathways_statistic2 <- aggregate(pathways$count, by = list(pathways$term), FUN = sum)
pathways_statistic <- data.frame(term = pathways_statistic$Group.1, 
                                na_count = pathways_statistic$x, 
                                total_count = pathways_statistic2$x) %>%
  mutate(na_fraction = na_count / total_count) %>%
  filter(na_fraction < 0.2) # remove pathways with over 20% genes failed to map

# filter pathways that lost too many genes during translation from HS to MM
pathways <- pathways %>%
  filter(term %in% pathways_statistic$term)

# convert pathway to list
pathway_signature_list <- list()
for (pathway in unique(pathways$term)) {
  data <- pathways %>% 
    filter(term == pathway) 
  pathway_signature_list[[pathway]] <- data$mm_gene
}

MB_integrated <- AddModuleScore(object = MB_integrated, 
                                 features = pathway_signature_list, 
                                 name = names(pathway_signature_list),
                                 assay = "RNA")

# !!!! remove the tailing numbers of pathways added by seurat !!!!
# need to modify if order of pathway name changed !!!!
colnames(MB_integrated@meta.data)[c(54:62)] <- gsub("\\d*$", "", colnames(MB_integrated@meta.data)[c(54:62)])

#####cell cluster annotation#####
#marker gene 
#https://www.nature.com/articles/s41467-019-13657-6
features_set <- c(
  "SHH-A", "SHH-B", "SHH-C",
  "Gfap", "Aqp4", "Slc1a3", #astrocyte
  "Atoh1", "Gli1", "Tubb3", "Sox11", # Tumor cell
  "Top2a", "Cdk1", "Rrm2", #SHH-A 
  "Eif3e", "Eef1a1", "Ptch1", "Boc",#SHH-B
  "Stmn2", "Map1b", "Tubb2b", "Sema6a", #SHH-C
  "Mki67", "Pcna", #cycling cell
  "Ptprc", # immune
  "Vwf", "Pecam1", "Cldn5",# endothelial
  "Vsnl1", "Scn1b", #neuron
  "Cd68", "Tyrobp", # myeloid  lineage
  "Pdgfrb", "Kcnj8", "Notch3", "Rgs5", #mural cell #https://www.nature.com/articles/nature25739
  "Mog", "Mag", # oligodendrocyte
  "Sox2","Pdgfra", "Olig1", "Olig2", "Sox10" #OPC
)

theme_my <- function (xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) 
{
  theme_classic() %+replace% theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = xlength, arrow = arrow), axis.line.y.left = element_line2(id = 2,
                                                                                                                                                         ylength = ylength, arrow = arrow), axis.title = element_text(hjust = 0.1))
}

dimplot_unannotated <- DimPlot(MB_integrated, label = T, reduction = "tsne", group.by = "integrated_snn_res.0.8", raster.dpi = c(1000, 1000), pt.size = 0.5) +
  scale_color_manual(values = colorRampPalette(brewer.pal(n = 10, name = "Paired"))(length(unique(MB_integrated$integrated_snn_res.0.8)))) + 
  theme_my()

featureplot_SHH_A <- FeaturePlot(MB_integrated, reduction = "tsne",features = c("SHH-A"),) +
  scale_colour_gradient(low = "#efefef", high = "red3") +
  theme_my() +
  theme(legend.key.size = unit(0.4, "cm"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6))

featureplot_SHH_B <- FeaturePlot(MB_integrated, reduction = "tsne",features = c("SHH-B"),) +
  scale_colour_gradient(low = "#efefef", high = "red3") +
  theme_my() +
  theme(legend.key.size = unit(0.4, "cm"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6))

featureplot_SHH_C <- FeaturePlot(MB_integrated, reduction = "tsne",features = c("SHH-C"),) +
  scale_colour_gradient(low = "#efefef", high = "red3") +
  theme_my() +
  theme(legend.key.size = unit(0.4, "cm"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6))

vlnplot_unannotated <- VlnPlot(MB_integrated, features = features_set, stack = T, flip = T, group.by = "integrated_snn_res.0.8", assay = "RNA") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(n = 10, name = "Paired"))(length(unique(features_set)))) + 
  NoLegend() + xlab("Cluster") +
  theme(strip.text.y = element_text(size = 9, face = "italic", hjust = 0),
        axis.text.x = element_text(size = 9, angle = 0))

dotplot_unannotated <- DotPlot(MB_integrated, features = rev(features_set), group.by = "integrated_snn_res.0.8", dot.scale = 6, assay = "RNA") +
  coord_flip() +  
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.3) +
  viridis::scale_colour_viridis(option="magma", direction = 1) +
  guides(size=guide_legend("Perc. exp.", override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.y = element_text(face = "italic", size = 7),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 0, size = 9))

p <- dimplot_unannotated + vlnplot_unannotated + dimplot_Sox2 +
  featureplot_SHH_A + featureplot_SHH_B + featureplot_SHH_C + plot_layout(nrow = 2, ncol = 3)
ggsave(plot = p, filename = "./plots/MB_integrated_cluster_dimplot_vlnplot_unanotated.pdf", width = 51, height = 30, units = "cm")

p <- dimplot_unannotated + dotplot_unannotated + dimplot_Sox2 +
  featureplot_SHH_A + featureplot_SHH_B + featureplot_SHH_C + plot_layout(nrow = 2, ncol = 3)
ggsave(plot = p, filename = "./plots/MB_integrated_cluster_dimplot_dotplot_unanotated.pdf", width = 51, height = 30, units = "cm")

# cluster identification
MB_integrated$cell_type <- plyr::mapvalues(MB_integrated$integrated_snn_res.0.8,
                                           from=0:22,
                                           to=c("Tumor cell",#0
                                                "Tumor cell",#1
                                                "Tumor cell",#2
                                                "Tumor cell",#3 
                                                "Tumor cell",#4
                                                "Tumor cell",#5
                                                "Tumor cell",#6
                                                "Tumor cell",#7
                                                "Neuron",#8
                                                "Tumor cell",#9
                                                "Tumor cell",#10
                                                "Tumor cell",#11
                                                "Tumor cell",#12
                                                "Astrocyte",#13
                                                "Myeloid cell",#14
                                                "Tumor cell",#15
                                                "Oligodendrocyte progenitor cell",#16
                                                "Neuron",#17
                                                "Neuron",#18
                                                "Pericytes",#19 
                                                "Oligodendrocyte",#20 
                                                "Oligodendrocyte",#21 
                                                "Endothelial cell"#22  
                                           ))
MB_integrated <- SetIdent(MB_integrated, value = MB_integrated$cell_type)
markers <- FindAllMarkers(MB_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "cell_type")
save(markers, file = "./processed_data/MB_integrated_annotated_cluster_markers.Rdata")
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
export::table2excel(top_markers, file = "./processed_data/MB_integrated_annotated_clusters_top_markers.xlsx")

# subset Pericyte cluster thet possiblely contain fibroblast
cluster_19 <- subset(MB_integrated, subset = (cell_type == "Pericytes"))
cluster_19 <- RunUMAP(cluster_19, reduction = "pca", dims = 1:14)
cluster_19 <- RunTSNE(cluster_19, reduction = "pca", dims = 1:14)
cluster_19 <- FindNeighbors(cluster_19, reduction = "pca", dims = 1:14)
DefaultAssay(cluster_19) <- "integrated"
cluster_19 <- FindClusters(cluster_19, resolution = seq(0.1,2,0.1))
clustree(cluster_19) +
  scale_color_manual(values = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(20)) +
  guides(color = guide_legend(ncol = 3)) #0.9 optimal
dimplot <- DimPlot(cluster_19, reduction = "tsne", group.by = "integrated_snn_res.0.9") +
  scale_color_manual(values = colorRampPalette(brewer.pal(n = 10, name = "Paired"))(length(unique(cluster_19$integrated_snn_res.0.9)))) +
  theme_my()

vlnplot <- VlnPlot(cluster_19, features = c("Pdgfra", "Vim", "Col1a1", "Col1a2", "Col5a1", "Loxl1", "Lum", "Fbln1", "Fbln2",  "Pdgfrb", "Rgs5", "Fbn1","Acta2", "Des", "Mcam", "Tagln", "Notch3", "Col3a1", "Cspg4", "Fap", "Pdpn", "S100a4", "Tnc", "Cd248"), stack = T, flip = T, group.by = "integrated_snn_res.0.8", assay = "RNA") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(n = 10, name = "Paired"))(24)) +
  NoLegend() + xlab("Cluster") +
  theme(strip.text.y = element_text(size = 6, face = "italic", hjust = 0),
        axis.text.x = element_text(size = 9, angle = 0))

dimplot + vlnplot
ggsave(plot = dimplot + vlnplot, filename = "./plots/QC/MB_integrated_cluster_19_pericyte.pdf", width = 20, height = 10, units = "cm")

# Cspg4 as pericyte marker
pericyte_id <- WhichCells(MB_integrated, idents = "Pericytes", slot = "count", expression = Cspg4 > 0)
fibroblast_id <- WhichCells(MB_integrated, idents = "Pericytes", slot = "count", expression = Cspg4 == 0)

cell_ident <- MB_integrated@meta.data[,c("cell_type", "orig.ident")] %>%
  mutate(cell_type = as.character(cell_type)) 
cell_ident$cell_type[rownames(cell_ident) %in% pericyte_id] <- "Pericyte"
cell_ident$cell_type[rownames(cell_ident) %in% fibroblast_id] <- "Fibroblast"
cell_ident$cell_type <- as.factor(cell_ident$cell_type)

MB_integrated$cell_type <- cell_ident$cell_type

#####Define SHH-A SHH-B SHH-C tumor cell#####
cell_type <- MB_integrated@meta.data[, c("SHH-A", "SHH-B", "SHH-C", "cell_type")]
cell_type$SHH_subtype <- NA

for (row in 1:nrow(cell_type)) {
  if (cell_type[row, "cell_type"] == "Tumor cell") {
    max <- max(cell_type[row, c("SHH-A", "SHH-B", "SHH-C")])
    cell_type$SHH_subtype[row] <- c("SHH-A", "SHH-B", "SHH-C")[which(cell_type[row, c("SHH-A", "SHH-B", "SHH-C")] == max)]
  }
}

table(cell_type$SHH_subtype)
# SHH-A SHH-B SHH-C 
# 3169 17601 29706

#add SHH subtype to object
MB_integrated$SHH_subtype <- cell_type$SHH_subtype
table(cell_type$SHH_subtype, MB_integrated$integrated_snn_res.0.8)

dimplot <- DimPlot(MB_integrated, reduction = "tsne", group.by = "SHH_subtype", pt.size = 0.01) +
  scale_colour_manual(values = c("#f3e721", "#2196f3", "#f3212d", "#efefef"), labels = c("SHH-A", "SHH-B", "SHH-C", "Other cell")) +
  theme_my() +
  theme(plot.title = element_blank(), 
        legend.key.size = unit(0.4, "cm"), 
        legend.title = element_blank(),  
        legend.text = element_markdown(size = 12)) # parsing markdown code for italic text
ggsave(plot = dimplot, filename = "./plots/MB_integrated_dimplot_SHH_subtype.pdf", width = 10, height = 7, units = "cm")

#####CytoTRACE#####
library(CytoTRACE)
Sys.setenv(RETICULATE_PYTHON="/data/Shared_software/anaconda3/envs/cytotrace/bin/python")

cytotrace_results <- CytoTRACE(as.matrix(MB_integrated@assays$RNA@counts), ncores = 128, subsamplesize = 1000)
cytotrace_results <- list(
  cytotrace_results[["CytoTRACE"]],
  cytotrace_results[["CytoTRACErank"]],
  cytotrace_results[["cytoGenes"]],
  cytotrace_results[["GCS"]],
  cytotrace_results[["gcsGenes"]],
  cytotrace_results[["Counts"]]
)
save(cytotrace_results, file = "./processed_data/MB_integrated_cytotrace_result.Rdata")

#add cytotrace res to object
MB_integrated$CytoTRACE <- cytotrace_results[[1]]
MB_integrated$CytoTRACErank <- cytotrace_results[[2]]

featureplot <- FeaturePlot(MB_integrated, features = c("CytoTRACE"), reduction = "tsne", pt.size = 0.01)  + 
  #scale_colour_gradientn(low = "#2d3694", mid = "lightyellow", high = "#ec342c", midpoint = 0.5) +
  scale_colour_gradientn(colors = c("#2f3690", "#2b85c9","#299672", "#c9c571", "#f69b2f", "#ee332b")) +
  theme_my() +
  theme(legend.key.size = unit(0.4, "cm"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6))
ggsave(plot = featureplot, filename = "./plots/MB_integrated_dimplot_cytotrace.pdf", width = 8, height = 7, units = "cm")

#save annotated object
save(MB_integrated, file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/MB_integrated.Rdata")

#####extract pathway score#####
pathway_score <- MB_integrated@meta.data %>%
  select(orig.ident, batch, genotype, sample_region, cell_type, S.Score, G2M.Score, Phase, Sox2_status, SHH_subtype, CytoTRACE, CytoTRACErank, seurat_clusters,
         names(pathway_signature_list),
         c("SHH-A", "SHH-B", "SHH-C")
         )
save(pathway_score, file = "./processed_data/MB_integrated_pathway_score.Rdata")

######plot SHH-A SHH-B SHH-C tumor cell proportion along differentiation path#####
library(ggridges)
library(patchwork)
library(ggpubr)
load("./processed_data/MB_integrated_pathway_score.Rdata")

# plot ridge plot
density_plot_by_region_Core_Skull <- pathway_score %>% 
  mutate(seurat_clusters = paste0("Cluster ",seurat_clusters)) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = SHH_subtype_cytotrace$seurat_clusters)) %>%
  filter(batch == "Experiment1") %>%
  filter(cell_type == "Tumor cell") %>%
  mutate(Region = ifelse(orig.ident == "WT_core_Experiment1", "Core tumor", "Peripheral tumor(Skull)")) %>%
  ggplot() +
  geom_density(aes(y = CytoTRACE, group = Region, fill = Region, color = Region), adjust = 1.5, alpha = .15) +
  scale_color_manual(values = c("#808080", "#0966cd")) + 
  scale_fill_manual(values = c("#808080", "#0966cd")) +
  labs(x = "Distribution density") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.text.x = element_text(size = 13))

# calculate statistics for box plot
stat.test_Core_Skull <- pathway_score %>% 
  mutate(seurat_clusters = paste0("Cluster ",seurat_clusters)) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = SHH_subtype_cytotrace$seurat_clusters)) %>%
  filter(batch == "Experiment1") %>%
  filter(cell_type == "Tumor cell") %>%
  mutate(Region = ifelse(orig.ident == "WT_core_Experiment1", "Core tumor", "Peripheral tumor(Skull)")) %>%
  t_test(CytoTRACE ~ Region) %>%
  add_xy_position() %>%  # add plot x y position
  mutate(p.signif = ifelse(p < 0.05, round(p, digits = 4), p.signif)) %>%
  mutate(p.signif = ifelse(p < 0.0001, "<0.0001", p.signif))

stat.test_Core_Skull$xmin <- -0.185 # fix alignment issue
stat.test_Core_Skull$xmax[which(stat.test_Core_Skull$xmax == 2)] <- 0.18
stat.test_Core_Skull$y.position <- stat.test_Core_Skull$y.position + 0.05

# boxplot
box_plot_by_region_Core_Skull <- pathway_score %>% 
  mutate(seurat_clusters = paste0("Cluster ",seurat_clusters)) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = SHH_subtype_cytotrace$seurat_clusters)) %>%
  filter(batch == "Experiment1") %>%
  filter(cell_type == "Tumor cell") %>%
  mutate(Region = ifelse(orig.ident == "WT_core_Experiment1", "Core tumor", "Peripheral tumor(Skull)")) %>%
  ggplot() +
  geom_boxplot(aes(y = CytoTRACE, group = Region, fill = Region),color = "black", outlier.shape = NA) +
  stat_pvalue_manual(stat.test_Core_Skull, label = "p.signif", size = 2.5) +
  scale_color_manual(values = c("#808080", "#0966cd")) + 
  scale_fill_manual(values = c("#808080", "#0966cd")) +
  labs(x = "") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.06))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.text.x = element_text(size = 13))

p <- density_plot_by_region_Core_Skull + box_plot_by_region_Core_Skull + 
  plot_layout(nrow = 2, ncol = 1)
ggsave(plot = p, filename = "./plots/MB_integrated_cytotrace_density_box_plot_by_region_Core_Skull.pdf", width = 3.4, height = 6)
