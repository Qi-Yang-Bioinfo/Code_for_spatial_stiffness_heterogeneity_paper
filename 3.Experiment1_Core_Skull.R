setwd("/home/oxygen/Ali_MB/")
library(infercnv)
library(Seurat) 
library(SeuratObject) 
library(patchwork)
library(ggtext)
library(tidyverse)
library(RColorBrewer)
library(tidydr)
library(data.table)
library(rstatix)
library(ggpubr)

source("./R_function/bezier_curve.R")

#####subset experiment1#####
# load("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/MB_integrated.Rdata")
# 
# Experiment1_Core_Skull <- subset(MB_integrated, subset = (orig.ident %in% c("WT_core_Experiment1", "WT_peripheral_adjacent_to_skull_Experiment1")))
# save(Experiment1_Core_Skull, file = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/MB_integrated_Experiment1_Core_Skull.Rdata")

load("/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/MB_integrated_Experiment1_Core_Skull.Rdata")

##### plot tsne dimension reduction #####
theme_my <- function (xlength = 0.25, ylength = 0.25, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) {theme_classic() %+replace% theme_noaxis(axis.line.x.bottom = element_line2(id = 1, xlength = xlength, arrow = arrow), axis.line.y.left = element_line2(id = 2, ylength = ylength, arrow = arrow), axis.title = element_text(hjust = 0.1))}

Experiment1_Core_Skull$Region <- ifelse(Experiment1_Core_Skull$orig.ident == "WT_core_Experiment1", "Core tumor", "Peripheral tumor(Skull)")

# find xlim ylim
tsne_x <- c(min(Experiment1_Core_Skull@reductions[["tsne"]]@cell.embeddings[,1]), max(Experiment1_Core_Skull@reductions[["tsne"]]@cell.embeddings[,1]))
tsne_y <- c(min(Experiment1_Core_Skull@reductions[["tsne"]]@cell.embeddings[,2]), max(Experiment1_Core_Skull@reductions[["tsne"]]@cell.embeddings[,2]))

region_dimplot <- DimPlot(Experiment1_Core_Skull, reduction = "tsne", group.by = "Region", pt.size = 0.05, shuffle = T,
                              cols = c("#acacac", "#0966cd") ) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 4), ncol=1)) + #smaller legend size
  theme(legend.text = element_text(size = 5), aspect.ratio = 1) +
  xlim(tsne_x) +
  ylim(tsne_y) +
  theme_void()
 
region_dimplot_tumor_cell <- DimPlot(subset(Experiment1_Core_Skull, subset = (cell_type %in% c("Tumor cell"))), reduction = "tsne", group.by = "Region", pt.size = 0.05, shuffle = T,
                          cols = c("#acacac", "#0966cd") ) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 4), ncol=1)) + #smaller legend size
  theme(legend.text = element_text(size = 5), aspect.ratio = 1) +
  xlim(tsne_x) +
  ylim(tsne_y) +
  theme_void()

cell_type_dimplot <- DimPlot(Experiment1_Core_Skull, reduction = "tsne", group.by = "cell_type", pt.size = 0.05, shuffle = T, 
                                   cols = brewer.pal(n = 12, name = "Paired")[1:9]) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 2.5), ncol=1)) + #smaller legend size
  theme(legend.position = "none", aspect.ratio = 1) +
  xlim(tsne_x) +
  ylim(tsne_y) +
  theme_void()

cell_type_dimplot_label <- DimPlot(Experiment1_Core_Skull, reduction = "tsne", group.by = "cell_type", pt.size = 0.05, label = T, label.size = 3, repel = T, shuffle = T,
                             cols = brewer.pal(n = 12, name = "Paired")[1:9]) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 2.5), ncol=1)) + #smaller legend size
  theme(legend.position = "none", aspect.ratio = 1) +
  xlim(tsne_x) +
  ylim(tsne_y) +
  theme_void()

cell_type_dimplot_legend <- DimPlot(Experiment1_Core_Skull, reduction = "tsne", group.by = "cell_type", pt.size = 0.05, shuffle = T,
                                   cols = brewer.pal(n = 12, name = "Paired")[c(9:10, 1:7)]) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 2.5), ncol=1)) + #smaller legend size
  theme(legend.text = element_text(size = 5), aspect.ratio = 1) +
  xlim(tsne_x) +
  ylim(tsne_y) +
  theme_void()

SHH_subtype_dimplot <- DimPlot(subset(Experiment1_Core_Skull, subset = (SHH_subtype %in% c("SHH-A", "SHH-B", "SHH-C"))), reduction = "tsne", group.by = "SHH_subtype", pt.size = 0.05, shuffle = T, 
                             cols = c("#f3e721", "#2196f3", "#f3212d")) +  
  theme_my(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  labs(title = "") + guides(color = guide_legend(override.aes = list(size = 2.5), ncol=1)) + #smaller legend size
  theme(legend.text = element_text(size = 5), aspect.ratio = 1) +
  xlim(tsne_x) +
  ylim(tsne_y) +
  theme_void()

cytotrace_featureplot <- FeaturePlot(Experiment1_Core_Skull, features = c("CytoTRACE"), reduction = "tsne", pt.size = 0.05)  + 
  scale_colour_gradientn(colors = c("#2f3690", "#2b85c9","#299672", "#c9c571", "#f69b2f", "#ee332b")) +
  theme_my() +
  theme(legend.key.size = unit(0.4, "cm"), 
        aspect.ratio = 1,
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6)) +
  xlim(tsne_x) +
  ylim(tsne_y) +
  theme_void()

plot <- region_dimplot + region_dimplot_tumor_cell + cell_type_dimplot + cell_type_dimplot_label +
  cell_type_dimplot_legend + SHH_subtype_dimplot + cytotrace_featureplot + plot_spacer() +
  plot_layout(nrow = 2, ncol = 4)

ggsave(plot = plot, filename = "./plots/MB_Experiment1_Core_Peri_Skull_dimplot.pdf", width = 23, height = 10)

##### Stem cell proliferation pathway score boxplot #####
options(scipen = 233)
load("./processed_data/MB_integrated_pathway_score.Rdata")

pathway_score <- pathway_score %>%
  filter(batch == "Experiment1") %>%
  filter(orig.ident %in%  c("WT_core_Experiment1", "WT_peripheral_adjacent_to_skull_Experiment1")) %>%
  mutate(Region = case_when(
    orig.ident == "WT_core_Experiment1" ~ "WT(core)",
    orig.ident == "WT_peripheral_adjacent_to_skull_Experiment1" ~ "WT(periphery)"
  )) 

stat.test <- pathway_score %>%
  filter(cell_type == "Tumor cell") %>%
  wilcox_test(GOBP_STEM_CELL_PROLIFERATION ~ Region) %>%
  add_xy_position() %>%  # add plot x y position
  add_significance() %>%
  mutate(p.signif = ifelse(p < 0.0001, "<0.0001", paste0("=",round(p, digits = 4)))) %>%
  mutate(label = ifelse(p.signif == "ns", "ns", paste0("italic(p)", p.signif))) 
stat.test$y.position = 2.85

pathway_score_vln_boxplot <- pathway_score %>%
  mutate(GOBP_STEM_CELL_PROLIFERATION = scale(GOBP_STEM_CELL_PROLIFERATION)) %>%
  filter(cell_type == "Tumor cell") %>%
  ggplot() +
  geom_violin(mapping = aes(x = Region, y = GOBP_STEM_CELL_PROLIFERATION, fill = Region), color = "white", width = 0.9, trim = F) +
  geom_boxplot(mapping = aes(x = Region, y = GOBP_STEM_CELL_PROLIFERATION, fill = Region), color = "white", width= 0.1, outlier.shape = NA, lwd = 0.2) +
  scale_fill_manual(values = c("#acacac", "#0966cd")) +
  geom_bracket(data = stat.test, type = "expression",
               tip.length = c(0.02, 0.02)
  ) +
  ylim(-2.3, 3) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        aspect.ratio = 2.4/1
  )
ggsave(plot = pathway_score_vln_boxplot, filename = "./plots/MB_Experiment1_Core_Peri_Skull_boxplot_GOBP_STEM_CELL_PROLIFERATION.pdf", width = 3.2, height = 2.8)