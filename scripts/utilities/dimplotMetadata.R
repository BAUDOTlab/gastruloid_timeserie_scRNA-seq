dimplotMetadata <- function(SO, metadata, title=NULL, reduction=NULL){
    
    if (metadata == "day") {
        use.colors <- colors.stage
        title <- paste0(SO@project.name, " : cells stages")
    } else if (length(grep("celltype", metadata)) != 0) {
        use.colors <- colors.celltype
        title <- if (is.null(title)) paste0(SO@project.name," : cell identities") else title
    } else {
        use.colors <- NULL
    }
    
    if (is.null(title)){
        title <- paste0(SO@project.name, " : ", metadata)
    }    
    
    if (is.null(reduction)){
        reduction <- grep("umap", names(SO@reductions), value = TRUE)[1]
    }
    
    Idents(SO) <- metadata
    
    plot <- DimPlot(SO,
                    pt.size = 1,
                    reduction = reduction,
                    label = T,
                    label.size = 6,
                    repel = T,
                    cols = use.colors[levels(Idents(SO))]) +
        NoLegend() +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text = element_blank(),
              # axis.ticks = element_blank()
              line = element_blank())

    print(plot)
}
