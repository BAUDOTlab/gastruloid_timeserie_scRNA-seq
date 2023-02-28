library(Seurat)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]


cat("\n\n-------- Start working on", dataset, "-----------------------\n\n")

top.pcs <- 30

cat("Load data\n")
dataset.anFolder <- file.path("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5", dataset)
SOsing <- readRDS(paste0("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/ciml5/rdsObjects/05_preprocessing_lab_4_days_merged.rds"))
atlas.subset <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/inputData/atlas/atlas_preprocessed.rds")

cat("\nPerform label transfer\n")
anchors <- FindTransferAnchors(reference = atlas.subset, query = SOsing,
                               dims = 1:top.pcs, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = atlas.subset$celltype,
                            dims = 1:top.pcs)

cellsData <- data.frame(predictions$predicted.id, SOsing$RNA_snn_res.1, row.names = rownames(SOsing@meta.data))
colnames(cellsData) <- c("multiTP_celltype", "multiTP_res1")

cat("\nWrite table\n")
write.table(
    x = cellsData,
    file = file.path(
        dataset.anFolder,
        paste0("metadata_celltypeANDclusters_", dataset, ".csv")),
    sep = ",",
    row.names = TRUE,
    col.names = TRUE)


# Rscript getMultipleMetadata.R lab_4_days
