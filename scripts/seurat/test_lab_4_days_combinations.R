library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)

SO4 <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/07_clusters_lab_day_04.rds")
SO5 <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/07_clusters_lab_day_05.rds")
SO6 <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/07_clusters_lab_day_06.rds")
SO11 <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/07_clusters_lab_day_11.rds")

invisible(sapply(c(SO4, SO5, SO6, SO11), function(SO){
    SO@meta.data <- SO@meta.data %>% select(-matches("^integrated_snn_res"))
    SO@meta.data <- SO@meta.data %>% select(-matches("^RNA_snn_res"))
    SO@meta.data <- SO@meta.data %>% select(-matches("^outlier."))
    SO@meta.data <- SO@meta.data %>% select(-matches("^DF.class"))
    SO@meta.data <- SO@meta.data %>% select(-matches("^pANN_"))
    SO@meta.data <- SO@meta.data %>% select(-matches("^seurat"))
}))

gc()

# Merge
SO <- merge(SO4, c(SO5, SO6, SO11))
SO <- ScaleData(SO, features=rownames(SO), do.scale=FALSE, verbose=FALSE)
SO <- FindVariableFeatures(SO, selection.method = "vst", nfeatures = 2000)
SO <- RunPCA(SO, npcs = 50, nfeatures.print = 5, seed.use = 17, verbose=TRUE)

saveRDS(SO, file = "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/00_lab_4_days_merged.rds")

p1 <- DimPlot(object = SO, reduction = "pca", pt.size = .1, group.by = "day")
p2 <- VlnPlot(object = SO, features = "PC_1", group.by = "day", pt.size = .1)
plot_grid(p1,p2)

# Harmony integration
harm <- RunHarmony(SO, "day", plot_convergence = TRUE)
harm <- ScaleData(harm, do.scale=FALSE, verbose = FALSE)
harm <- RunPCA(harm, npcs = 50, verbose = FALSE)

saveRDS(harm, file = "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/00_lab_4_days_harmony.rds")

p1 <- DimPlot(object = harm, reduction = "harmony", pt.size = .1, group.by = "day")
p2 <- VlnPlot(object = harm, features = "harmony_1", group.by = "day", pt.size = .1)
plot_grid(p1,p2)

rm(harm, SO)
gc()

# Seurat block integration
SO.list <- list(SO4, SO5, SO6, SO11)
features <- SelectIntegrationFeatures(object.list = SO.list)
anchors <- FindIntegrationAnchors(object.list = SO.list, anchor.features = features)
blkS <- IntegrateData(anchorset = anchors)

blkS <- ScaleData(blkS, do.scale=FALSE, verbose = FALSE)
blkS <- RunPCA(blkS, npcs = 50, verbose = FALSE)

saveRDS(blkS, file = "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/00_lab_4_days_seuratBlk.rds")

p1 <- DimPlot(object = blkS, reduction = "pca", pt.size = .1, group.by = "day")
p2 <- VlnPlot(object = blkS, features = "PC_1", group.by = "day", pt.size = .1)
plot_grid(p1,p2)

# Seurat sequential integration
seqS <- SO.list[[1]]
names(SO.list) <- c("Days04", "Day_05", "Day_06", "Day_11")
for (a_day in names(SO.list)[2:length(SO.list)]){
    integ.list <- list(seqS, SO.list[[a_day]])
    features <- SelectIntegrationFeatures(integ.list)
    anchors <- FindIntegrationAnchors(object.list = integ.list,
                                      anchor.features = features)
    seqS <- IntegrateData(anchorset = anchors)
} 

seqS <- ScaleData(seqS, do.scale=FALSE, verbose = FALSE)
seqS <- RunPCA(seqS, npcs = 50, verbose = FALSE)

p1 <- DimPlot(object = seqS, reduction = "pca", pt.size = .1, group.by = "day")
p2 <- VlnPlot(object = seqS, features = "PC_1", group.by = "day", pt.size = .1)
plot_grid(p1,p2)

saveRDS(seqS, file = "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/00_lab_4_days_seuratSeq.rds")
