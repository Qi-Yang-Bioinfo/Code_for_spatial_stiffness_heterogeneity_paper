setwd("/home/oxygen/Ali_MB/")
library(Seurat) 
library(SeuratObject) 
library(patchwork)
library(ggtext)
library(tidyverse)
library(RColorBrewer)
library(tidydr)
library(data.table)
library(ComplexHeatmap)
library(circlize)

##### plot cell proportion Experiment1 #####
load("./processed_data/MB_integrated_pathway_score.Rdata")

cell_count <- pathway_score %>% 
  filter(batch == "Experiment2") %>%
  filter(orig.ident %in%  c("WT_core_Experiment2", "WT_peripheral_Experiment2", "Gpr68KO_core_Experiment2", "Gpr68KO_peripheral_Experiment2")) %>%
  mutate(Region = case_when(
    orig.ident == "WT_core_Experiment2" ~ "WT(core)",
    orig.ident == "WT_peripheral_Experiment2" ~ "WT(periphery)",
    orig.ident == "Gpr68KO_core_Experiment2" ~ "Gpr68KO(core)",
    orig.ident == "Gpr68KO_peripheral_Experiment2" ~ "Gpr68KO(periphery)"
  )) %>%
  select(orig.ident, Region, seurat_clusters, SHH_subtype, cell_type) 
  
cell_count_propotion <- data.frame(prop.table(table(cell_count$Region, cell_count$cell_type), margin = 1))
colnames(cell_count_propotion) = c("Region", "cell_type", "value")
cell_count_propotion_plot <- ggplot(cell_count_propotion, aes(Region, value, fill = cell_type)) +
  geom_col() + xlab("") + ylab("Proportion of cells (%)") +
  scale_fill_manual(values = brewer.pal(n = 12, name = "Paired")[c(9:10, 1:6)]) + 
  labs(fill = "Cell type") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 0.5, b = 0, l = 10, unit = "pt"),
    aspect.ratio = 3/2
  )

cell_count <- data.frame(table(cell_count$Region, cell_count$cell_type))
colnames(cell_count) = c("Region", "cell_type", "value")
cell_count_plot <- ggplot(cell_count, aes(Region, value, fill = cell_type)) +
  geom_col() + xlab("") + ylab("Cell number") +
  scale_fill_manual(values = brewer.pal(n = 12, name = "Paired")[c(9:10, 1:6)]) + 
  labs(fill = "Cell type") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(),
        plot.margin = margin(t = 5, r = 0, b = 0, l = 5, unit = "pt"),
        aspect.ratio = 3/2)

p <- cell_count_propotion_plot + cell_count_plot

ggsave(plot = p, filename = "./plots/MB_Experiment2_stacked_barplot_count_of_cell_type.pdf", height = 3.1, width = 10, units = "in")
export::table2excel(cbind(cell_count, cell_count_propotion), file = "./processed_data/MB_Experiment2_cell_count_prop.xlsx")

##### Stress fiber, Focal adhesion, Response to mechanical stimulus, Cilia, SHH signaling, Proliferation, Stemness #####
#pathway score heatmao
options(scipen = 233)
load("./processed_data/MB_integrated_pathway_score.Rdata")

# list pathways
selected_pathway <- c(
  
  #Mechanic stimulus response
  "GOBP_RESPONSE_TO_MECHANICAL_STIMULUS", 
  "GOBP_CELLULAR_RESPONSE_TO_MECHANICAL_STIMULUS", 
  
  #SHH pathway
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
    
  #Stress fiber
  "GOBP_STRESS_FIBER_ASSEMBLY", 
  "GOBP_POSITIVE_REGULATION_OF_STRESS_FIBER_ASSEMBLY", 
    
  #proliferation
  "GOBP_STEM_CELL_PROLIFERATION",
  "GOBP_REGULATION_OF_STEM_CELL_PROLIFERATION",
  "Stemness"
)

# filter pathway
pathway_score <- pathway_score %>%
  filter(orig.ident %in%  c("WT_core_Experiment1", "WT_peripheral_adjacent_to_skull_Experiment1", "Gpr68KO_core_Experiment2", "Gpr68KO_peripheral_Experiment2")) %>%
  mutate(Region = case_when(
    orig.ident == "WT_core_Experiment1" ~ "WT core",
    orig.ident == "WT_peripheral_adjacent_to_skull_Experiment1" ~ "WT periphery",
    orig.ident == "Gpr68KO_core_Experiment2" ~ "Gpr68KO core",
    orig.ident == "Gpr68KO_peripheral_Experiment2" ~ "Gpr68KO periphery"
  )) %>%
  filter(cell_type == "Tumor cell") %>%
  select(all_of(selected_pathway), Region) 

pathway_score_mean <- aggregate(pathway_score[,c(selected_pathway)], by = list(pathway_score$Region), FUN = mean) %>%
  column_to_rownames("Group.1") %>%
  scale() %>%
  t() %>% as.data.frame()

pathway_score_mean <- pathway_score_mean[, c("WT core", "WT periphery", "Gpr68KO core", "Gpr68KO periphery")]

colnames(pathway_score_mean) <- c("Core", "Periphery", "Core", "Periphery")

rownames(pathway_score_mean) <-  c(
  #Mechanic stimulus response
  "GOBP Response to mechanical stimulus", 
  "GOBP Cellular response to mechanical stimulus", 
  
  #SHH pathway
  "KEGG Hedgehog signaling pathway",
  
  #Stress fiber
  "GOBP Stress fiber assembly", 
  "GOBP Positive regulation of stress fiber assembly", 
  
  #proliferation
  "GOBP Stem cell proliferation",
  "GOBP Regulation of stem cell proliferation",
  "Stemness signature (Alex et al. PNAS 2019)"
)

pathway_split <- c(
  rep("Mechano-response", 2),
  rep("Hedgehog signaling", 1),
  rep("Stress fiber", 2),
  rep("Stem cell proliferation", 3)
)

pathway_split <- factor(x = pathway_split, levels = c("Mechano-response", "Hedgehog signaling", "Stress fiber", "Focal adhesion", "Stem cell proliferation"))

max <- max(pathway_score_mean)
min <- min(pathway_score_mean)
middle <- (max + min) / 2

heatmap <- Heatmap(pathway_score_mean,
 column_split = c("WT", "WT", "Gpr68 KO", "Gpr68 KO"),
 row_split = pathway_split,
 row_title = NULL,
 show_column_dend = F,
 show_column_names = T,
 show_row_names = T,
 column_names_rot = 45,
 cluster_rows = F,
 cluster_columns = F,
 width = unit(3.5, "cm"), 
 height = unit(7, "cm"),
 name = "Mean score",
 cell_fun = function(j, i, x, y, width, height, fill) {
 grid.text(sprintf("%.2f", pathway_score_mean[i, j]), x, y, gp = gpar(fontsize = 8))
 },
 rect_gp = gpar(col = "black", lwd = 1),
 heatmap_legend_param = list(border = "black", direction = "horizontal", legend_height = unit(1.5, "cm")),
 col = colorRamp2(c(min, middle, max), c("#4c588a", "white", "#de2e29")))

heatmap_no_text <- Heatmap(pathway_score_mean,
                           column_split = c("WT", "WT", "Gpr68 KO", "Gpr68 KO"),
                           row_split = pathway_split,
                   row_title = NULL,
                   show_column_dend = F,
                   show_column_names = T,
                   column_names_rot = 45,
                   show_row_names = T,
                   cluster_rows = F,
                   cluster_columns = F,
                   width = unit(3.5, "cm"), 
                   height = unit(7, "cm"),
                   name = "Mean score",
                   rect_gp = gpar(col = "black", lwd = 1),
                   heatmap_legend_param = list(border = "black", direction = "horizontal", legend_height = unit(1.5, "cm")),
                   col = colorRamp2(c(min, middle, max), c("#4c588a", "white", "#de2e29")))

pdf(file = "./plots/MB_Experiment1_2_pathway_score_heatmap.pdf", width = 12, height = 6)
draw(heatmap + heatmap_no_text, heatmap_legend_side = "top")
dev.off()
