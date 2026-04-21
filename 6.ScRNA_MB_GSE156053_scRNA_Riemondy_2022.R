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

#####preprocess#####
# load data
mat <- data.table::fread("/data/Sequencing_data_inhouse/Medulloblastoma/GSE156053_Human_MB_scRNA_Riemondy_2022/GSE155446/GSE155446_human_raw_counts.csv.gz")
mat <- data.frame(mat[,-1], row.names = gsub(".+[|]", "", mat[,1][[1]]))

meta <- data.table::fread("/data/Sequencing_data_inhouse/Medulloblastoma/GSE156053_Human_MB_scRNA_Riemondy_2022/GSE155446/GSE155446_human_cell_metadata.csv.gz") %>%
  as.data.frame() %>%
  column_to_rownames("cell")

#create seurat object
GSE156053_Human_MB_scRNA_Riemondy_2022 <- CreateSeuratObject(counts = mat,
                               project = "GSE156053_Human_MB_scRNA_Riemondy_2022",
                               meta.data = meta)

# preprocess
GSE156053_Human_MB_scRNA_Riemondy_2022 <- NormalizeData(GSE156053_Human_MB_scRNA_Riemondy_2022, normalization.method = "LogNormalize", scale.factor = 10000)
GSE156053_Human_MB_scRNA_Riemondy_2022 <- FindVariableFeatures(GSE156053_Human_MB_scRNA_Riemondy_2022, selection.method = "vst", nfeatures = 2000)
GSE156053_Human_MB_scRNA_Riemondy_2022 <- ScaleData(GSE156053_Human_MB_scRNA_Riemondy_2022)
GSE156053_Human_MB_scRNA_Riemondy_2022 <- RunPCA(GSE156053_Human_MB_scRNA_Riemondy_2022, npcs = 30, verbose = FALSE)
Elbow_plots <- ElbowPlot(GSE156053_Human_MB_scRNA_Riemondy_2022) + labs(title = GSE156053_Human_MB_scRNA_Riemondy_2022@project.name) &
  theme(
    plot.title = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),
    axis.title = element_text(size = 6),
    axis.line = element_line(size = 0.1),
    axis.ticks = element_line(size = 0.1)
  )
export::graph2pdf(x = Elbow_plots, file = paste0("./plots/QC/GSE156053_Human_MB_scRNA_Riemondy_2022_Elbow_plots.pdf"), width = 2, height = 2)
GSE156053_Human_MB_scRNA_Riemondy_2022 <- FindNeighbors(GSE156053_Human_MB_scRNA_Riemondy_2022, dims = 1:10)
GSE156053_Human_MB_scRNA_Riemondy_2022 <- FindClusters(GSE156053_Human_MB_scRNA_Riemondy_2022, resolution = 10)
GSE156053_Human_MB_scRNA_Riemondy_2022 <- RunUMAP(GSE156053_Human_MB_scRNA_Riemondy_2022, dims = 1:10, verbose = FALSE)

theme_my <- function(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) { theme_classic() %+replace% theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = xlength, arrow = arrow), axis.line.y.left = element_line2(id = 2, ylength = ylength, arrow = arrow), axis.title = element_text(hjust = 0.1))}
dimplot <- DimPlot(GSE156053_Human_MB_scRNA_Riemondy_2022, reduction = "umap", group = "coarse_cell_type",
                   cols = colorRampPalette(brewer.pal(n = 11, name = "Paired"))(length(unique(GSE156053_Human_MB_scRNA_Riemondy_2022$coarse_cell_type))) )+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "") +
  theme_my()

dimplot_subgroup <- DimPlot(GSE156053_Human_MB_scRNA_Riemondy_2022, reduction = "umap", group = "subgroup",
                            cols = colorRampPalette(brewer.pal(n = 11, name = "Paired"))(length(unique(GSE156053_Human_MB_scRNA_Riemondy_2022$subgroup))) )+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "") +
  theme_my()

p <- dimplot + dimplot_subgroup
ggsave(filename = "./plots/GSE156053_Human_MB_scRNA_Riemondy2022_dimplot.pdf", plot = p, width = 12, height = 4.5)

# save preprocessed object
save(GSE156053_Human_MB_scRNA_Riemondy_2022, file = "/data/Sequencing_data_inhouse/Medulloblastoma/GSE156053_Human_MB_scRNA_Riemondy_2022/GSE156053_Human_MB_scRNA_Riemondy_2022_preprocessed.Rdata")

# define mechanosensors
Mechanosensor <- c("GPR68")

# plot expr of mechanosenser
theme_my <- function(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) { theme_classic() %+replace% theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = xlength, arrow = arrow), axis.line.y.left = element_line2(id = 2, ylength = ylength, arrow = arrow), axis.title = element_text(hjust = 0.1))}

plot_list <- lapply(Mechanosensor, function (x) {
  
  expr <- FetchData(GSE156053_Human_MB_scRNA_Riemondy_2022, vars = paste0(x), layer = "count")
  GSE156053_Human_MB_scRNA_Riemondy_2022$expr_positive_negative <- ifelse(expr > 0, "Positive", "Negative")
  
  DimPlot(GSE156053_Human_MB_scRNA_Riemondy_2022, group.by = c("expr_positive_negative")) +
    scale_colour_manual(values = c("#efefef", "red3")) +
    labs(title = x) +
    theme_my() +
    theme(plot.title = element_text(face = "italic"), 
          legend.key.size = unit(0.4, "cm"), 
          legend.title = element_text(size = 7), 
          legend.text = element_text(size = 6))
})

design <- "A#
           BC"

p <- plot_list[[1]] + 
  dimplot + dimplot_subgroup +
  plot_layout(design = design)

ggsave(filename = "./plots/GSE156053_Human_MB_scRNA_Riemondy2022_dimplot_mechanosenser_expr.pdf", plot = p, width = 6, height = 6)



