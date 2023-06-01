library(Seurat)
library(ggplot2)

source("../utilities/seuratToURD.R")

SO <- readRDS("../../../seuratAnalysis/ciml5/rdsObjects/99_rawData_filtered_lab_4_days.rds")

umap2d <- read.table("../../../seuratAnalysis/ciml5/lab_4_days_merged/final_umap2Dcoordinates_lab_4_days_merged.csv", header = TRUE, sep = ",", row.names = 1)
SO[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap2d), key = "UMAP_")

metadata <- read.table("../../../seuratAnalysis/ciml5/lab_4_days_merged/metadata_celltypeANDclusters_lab_4_days_merged.csv", header = TRUE, sep = ",", row.names = 1)
md04 <- read.table("../../../seuratAnalysis/ciml5/lab_day_04/metadata_celltypeANDclusters_lab_day_04.csv", header = TRUE, sep = ",", row.names = 1)
md05 <- read.table("../../../seuratAnalysis/ciml5/lab_day_05/metadata_celltypeANDclusters_lab_day_05.csv", header = TRUE, sep = ",", row.names = 1)
md06 <- read.table("../../../seuratAnalysis/ciml5/lab_day_06/metadata_celltypeANDclusters_lab_day_06.csv", header = TRUE, sep = ",", row.names = 1)
md11 <- read.table("../../../seuratAnalysis/ciml5/lab_day_11/metadata_celltypeANDclusters_lab_day_11.csv", header = TRUE, sep = ",", row.names = 1)

md04['day'] <- 'Day_04'
md05['day'] <- 'Day_05'
md06['day'] <- 'Day_06'
md11['day'] <- 'Day_11'


l <- list(md04, md05, md06, md11)
df1 <- Reduce(rbind, l)
df1 <- data.frame(df1, rn = row.names(df1))
SO@meta.data <- data.frame(SO@meta.data, rn = row.names(SO@meta.data))
df2 <- merge(SO@meta.data, df1, by = "rn", all = TRUE)
metadata2  <- data.frame(metadata, rn = row.names(metadata))
df3 <- merge(df2, metadata2, by = "rn", all = TRUE)
df4 <- transform(df3, row.names = rn, rn = NULL)

SO@meta.data <- df4

# Difficulties to install URD in this conda environment
# I simply copied the format of an URD object here, from urd-class.R
URD <- methods::setClass("URD", slots = c(
    count.data=c("dgCMatrix", NULL), 
    logupx.data=c("dgCMatrix", NULL), 
    meta="data.frame", 
    group.ids="data.frame", 
    var.genes="vector", 
    knn="list",
    pca.sdev="vector", 
    pca.load="data.frame", 
    pca.scores="data.frame", 
    pca.sig="vector", 
    tsne.y="data.frame", 
    plot.3d="list",
    gene.sig.z="data.frame", 
    dm=c("DiffusionMap",NULL), 
    diff.data="data.frame", 
    pseudotime="data.frame",
    pseudotime.stability="list", 
    tree="list",
    nmf.g=c("dgCMatrix", "matrix", NULL),
    nmf.c=c("dgCMatrix", "matrix", NULL)
))


urdO <- seuratToURD(SO)
saveRDS(urdO, file = "../../../lineageInference/inputData/lab_data/urdObject_lab_4_days.rds")
