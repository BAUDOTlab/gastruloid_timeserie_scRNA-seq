# Modified to fit the Seurat V3 format
# The seurat object has changed its format, so modifications make this function
# able to take data from a more recent Seurat object (v3 and v4).

#' Import URD from Seurat
#' 
#' If you're previously performed analyses in Seurat, you can copy them over
#' directly into an URD object for analyzing pseudotime, pseudolineage, and
#' building a specification tree.
#' 
#' @param seurat.object A Seurat object
#' @return An URD object
#' @export

seuratToURD <- function(seurat.object) {
    if (requireNamespace("Seurat", quietly = TRUE)) {
        # Create an empty URD object
        ds <- new("URD")
        
        # Copy over data
        ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
        if(!any(dim(seurat.object@assays$RNA@data) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$RNA@data[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
        
        # Copy over metadata
        ## TO DO - grab kmeans clustering info
        get.data <- NULL
        if (.hasSlot(seurat.object, "data.info")) { 
            get.data <- as.data.frame(seurat.object@assays$RNA@data.info)
        } else if (.hasSlot(seurat.object, "meta.data")) { 
            get.data <- as.data.frame(seurat.object@meta.data) 
        }
        if(!is.null(get.data)) {
            di <- colnames(get.data)
            m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
            discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
            gi <- di[which(discrete <= 0.015)]
            ds@meta <- get.data[,m,drop=F]
            ds@group.ids <- get.data[,gi,drop=F]
        }
        
        # Copy over var.genes
        if(length(seurat.object@assays$RNA@var.features > 0)) ds@var.genes <- seurat.object@assays$RNA@var.features
        
        # Move over tSNE projection
        if (.hasSlot(seurat.object, "tsne.rot")) {
            if(!any(dim(seurat.object@tsne.rot) == 0)) {
                ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
                colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
            }
        } else if (.hasSlot(seurat.object, "reductions")) {
            if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
                ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
                colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
            }
        }
        
        # Move over PCA results
        if (.hasSlot(seurat.object, "pca.x")) {
            if(!any(dim(seurat.object@pca.x) == 0)) {
                ds@pca.load <- seurat.object@pca.x
                ds@pca.scores <- seurat.object@pca.rot
                warning("Need to set which PCs are significant in @pca.sig")
            }
            ## TO DO: Convert SVD to sdev
        } else if (.hasSlot(seurat.object, "reductions")) {
            if(("pca" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$pca@gene.loadings) == 0)) {
                ds@pca.load <- as.data.frame(seurat.object@reductions$pca@gene.loadings)
                ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
                ds@pca.sdev <- seurat.object@reductions$pca@sdev
                ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
            }
        }
        return(ds)
    } else {
        stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
    }
}
