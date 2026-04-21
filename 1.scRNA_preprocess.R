setwd("/home/oxygen/Ali_MB/")

library(Seurat)
library(SeuratObject) 
library(patchwork)
library(harmony)
library(tidyverse)
library(DoubletFinder)
source("./R_function/ScRNAseq_preprocess_function.R")

#####Preprocess#####
# Inhouse data 
WT_core_Experiment1 <- readH5_QC_filter_doublet_remove(h5_mat_path = "/data/Sequencing_data_inhouse/Huang_lab/20231231_Ali_MB_snRNAseq/processed_data/Experiment1/WT_core/outs/filtered_feature_bc_matrix.h5",
                                        sample_id = "WT_core_Experiment1", doublet_rate = 0.08)
save(WT_core_Experiment1, file = paste0("/data/Sequencing_data_inhouse/Huang_lab/20231231_Ali_MB_snRNAseq/processed_data/Experiment1/WT_core/preprocessed_doublet_filtered_WT_core_Experiment1.Rdata"))

WT_peripheral_adjacent_to_skull_Experiment1 <- readH5_QC_filter_doublet_remove(h5_mat_path = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Experiment1/WT_peripheral_adjacent_to_skull/outs/filtered_feature_bc_matrix.h5",
                                               sample_id = "WT_peripheral_adjacent_to_skull_Experiment1", doublet_rate = 0.056)
save(WT_peripheral_adjacent_to_skull_Experiment1, file = paste0("/data/Sequencing_data_inhouse/Huang_lab/20231231_Ali_MB_snRNAseq/processed_data/Experiment1/WT_peripheral_adjacent_to_skull/preprocessed_doublet_filtered_WT_peripheral_adjacent_to_skull_Experiment1.Rdata"))

Gpr68KO_core_Experiment2 <- readH5_QC_filter_doublet_remove(h5_mat_path = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Experiment2/Gpr68KO_core/outs/filtered_feature_bc_matrix.h5",
                                   sample_id = "Gpr68KO_core_Experiment2", doublet_rate = 0.056)
save(Gpr68KO_core_Experiment2, file = paste0("/data/Sequencing_data_inhouse/Huang_lab/20231231_Ali_MB_snRNAseq/processed_data/Experiment2/Gpr68KO_core/preprocessed_doublet_filtered_Gpr68KO_core_Experiment2.Rdata"))

Gpr68KO_peripheral_Experiment2 <- readH5_QC_filter_doublet_remove(h5_mat_path = "/data/Sequencing_data_inhouse/Huang_lab/Ali_MB_snRNAseq/processed_data/Experiment2/Gpr68KO_peripheral/outs/filtered_feature_bc_matrix.h5",
                                   sample_id = "Gpr68KO_peripheral_Experiment2", doublet_rate = 0.056)
save(Gpr68KO_peripheral_Experiment2, file = paste0("/data/Sequencing_data_inhouse/Huang_lab/20231231_Ali_MB_snRNAseq/processed_data/Experiment2/Gpr68KO_peripheral/preprocessed_doublet_filtered_Gpr68KO_peripheral_Experiment2.Rdata"))