readH5_QC_filter_doublet_remove <-
  function(h5_mat_path = NULL, gz_mat_path = NULL, sample_id = NULL, doublet_rate = NULL) {
    
    # create object, plot QC 
    if (!is.null(h5_mat_path) & is.null(gz_mat_path) ) {
      print("Read H5 expr matrix!")
      data <- CreateSeuratObject(counts = Read10X_h5(h5_mat_path), project = sample_id, min.cells = 0, min.features = 0)
    }
    
    if (!is.null(gz_mat_path) & is.null(h5_mat_path)) {
      print("Read gziped expr matrix!")
      data <- CreateSeuratObject(counts = Read10X(gz_mat_path), project = sample_id, min.cells = 0, min.features = 0)
    }
    
    if (!is.null(h5_mat_path) & !is.null(gz_mat_path) ) {
      print("H5 and gzipped mat found, use H5 matrix!")
      data <- CreateSeuratObject(counts = Read10X_h5(h5_mat_path), project = sample_id, min.cells = 0, min.features = 0)
    }
    
    if (is.null(gz_mat_path) & is.null(h5_mat_path)) {
      stop("No matrix found!")
    }

    # add QC
    data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt") #mitochonrial
    data <- PercentageFeatureSet(data, pattern = "^RP[SL]", col.name = "percent.rb") #ribosome
    data <- PercentageFeatureSet(data, features = CaseMatch(c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"), rownames(data)), col.name = "percent.HB")   #erythrocyte
    
    Pre_filter_QC_plots <-
      VlnPlot(
        data,
        c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB"),
        pt.size = 0.001,
        stack = T,
        flip = T
      ) &
      theme(
        text = element_text(size = 6),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        legend.position = "none")
    
    #export::graph2office(x = Pre_filter_QC_plots,  file = paste0("./plots/QC/", sample_id,"_Pre_filter_QC_plots.pptx"), width = 1.5, height = 2)
    export::graph2pdf(x = Pre_filter_QC_plots, file = paste0("./plots/QC/", sample_id,"_Pre_filter_QC_plots.pdf"), width = 1.5, height = 2)
    export::graph2png(x = Pre_filter_QC_plots,  file = paste0("./plots/QC/", sample_id,"_Pre_filter_QC_plots.png"), width = 1.5, height = 2)
                      
    #filter by nFeature_RNA nCount_RNA percent.mt percent.rb
    data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA <7500 & percent.mt < 25 & percent.rb < 40 & percent.HB < 1)
    
    # preprocess
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures =2000)
    data <- ScaleData(data)
    data <- RunPCA(data, npcs = 30, verbose = FALSE)
    data <- FindNeighbors(data, dims = 1:15)
    data <- FindClusters(data, resolution = 0.5)
    data <- RunUMAP(data, dims = 1:15, verbose = FALSE)
    Elbow_plots <- ElbowPlot(data) + labs(title = data@project.name) &
      theme(
        plot.title = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1)
      )
    #export::graph2office(x = Elbow_plots, file = paste0("./plots/QC/", sample_id,"Elbow_plots.pptx"), width = 2, height = 2)
    export::graph2pdf(x = Elbow_plots, file = paste0("./plots/QC/", sample_id,"_Elbow_plots.pdf"), width = 2, height = 2)
    export::graph2png(x = Elbow_plots, file = paste0("./plots/QC/", sample_id,"_Elbow_plots.png"), width = 2, height = 2)
    
    # filter doublet
    print(paste0("DoubletFinder processing sample: ", data@project.name))
    #pK Identification (no ground-truth)
    sweep.res <- paramSweep(data, PCs = 1:15, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pk <- as.numeric(as.vector(bcmvn$pK)[which.max(bcmvn$BCmetric)])
    #Homotypic Doublet Proportion Estimate
    homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)  
    #calculate doublet rate
    doublet_rate <- doublet_rate 
    nExp_poi <- round(doublet_rate*nrow(data@meta.data))  
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    ## Run DoubletFinder with varying classification stringencies
    data <- doubletFinder(data, PCs = 1:15, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    data <- doubletFinder(data, PCs = 1:15, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_", pk, "_", nExp_poi), sct = FALSE)
    #rename metadata colname
    tmp_object_meta <- data@meta.data
    colnames(tmp_object_meta)[ncol(tmp_object_meta) - 1] <- "DoubletFinder_res_initial_run"
    colnames(tmp_object_meta)[ncol(tmp_object_meta)] <- "DoubletFinder_res_final"
    tmp_object_meta -> data@meta.data
    # remove doublet
    data <- subset(data, DoubletFinder_res_final == "Singlet")
    # clustering
    data <- FindNeighbors(data, reduction = "pca", dims = 1:15)
    data <- FindClusters(data, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

    
    # return object
    return(data)
  }
