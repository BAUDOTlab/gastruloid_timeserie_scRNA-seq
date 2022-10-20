library(Seurat)
library(SeuratObject)
library(SeuratDisk)

myargs <- commandArgs(trailingOnly=TRUE)

inputFile <- myargs[1] # ends by .rds
outputFile <- myargs[2] # ends with .h5Seurat

cat("\n-----------------------------------------------------------------------------\n")
print(paste0("Conversion of the file ", inputFile))
cat("\n-----------------------------------------------------------------------------\n")

#' Load the Seurat Object
cat("Load Seurat Object (SO)\n\n")
SO <- readRDS(inputFile) 

#' update the Seurat Object to at least v3.1.2
#cat("Update the SO\n\n")
#SO_up <- SeuratObject::UpdateSeuratObject(SO)

#' Save as a h5Seurat file format
cat("Save SO as a h5Seurat file\n\n")
SeuratDisk::SaveH5Seurat(SO, outputFile)

#' Save update RDS file
#fileUp <- SeuratDisk::LoadH5Seurat(outputFile)
#saveRDS(fileUp, file = paste0(outputFile, ".rds"))

#' #' Check the content of the H5Seurat file
#' SO_H5 <- SeuratDisk::LoadH5Seurat("/mnt/DATA_4TB/projects/gastruloids_scRNAseq/CardiacGastruloidsAnalysis/gatherAnalysis/rdsObjects/08_cellIdentity2LouvainClusters_zaffran3cumulates_alone.h5Seurat")

#' Convert to a h5ad file format, readable by scanpy
cat("Do the conversion from h5Seurat to h5ad file format\n\n")
SeuratDisk::Convert(outputFile, dest = "h5ad", assay = "RNA", overwrite = TRUE)
