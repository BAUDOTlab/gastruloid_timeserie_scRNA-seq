#!/usr/bin/env Rscript

dotPlotsMarkers <- function(SO, geneList, title=NULL){
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    
    ## DotPlot the n-top markers of each cluster
    nbMarkers <- length(geneList)
    if (nbMarkers %/% 40 == 0) {
        DotPlot(SO, features = geneList) +
            RotatedAxis() +
            CenterTitle() +
            ggtitle(paste0(title, "\n", SO@project.name)) +
            scale_colour_gradientn(colours = myPalette(100))
    } else {
        nbPlot <- (nbMarkers%/%40)+1
        for (i in seq(1, nbPlot)){
            print(DotPlot(SO, features = geneList[seq(1+(i-1)*(ceiling(nbMarkers/nbPlot)), i*ceiling(nbMarkers/nbPlot))]) +
                      RotatedAxis() +
                      CenterTitle() +
                      ggtitle(paste0(title, " - ", i, "/", (nbMarkers%/%40)+1, "\n", SO@project.name)) +
                      scale_colour_gradientn(colours = myPalette(100)))
        }
    }
}


dotPlotsGoi <- function(SO, geneList, title=NULL){
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    
    ## DotPlot the genes of interest
    nbMarkers <- length(geneList)
    if (nbMarkers %/% 40 == 0) {
        DotPlot(SO, features = sort(geneList)) +
            RotatedAxis() +
            CenterTitle() +
            ggtitle(paste0(title, "\n", SO@project.name)) +
            scale_colour_gradientn(colours = myPalette(100))
    } else {
        nbPlot <- (nbMarkers%/%40)+1
        for (i in seq(1, nbPlot)){
            print(DotPlot(SO, features = geneList[seq(1+(i-1)*(ceiling(nbMarkers/nbPlot)), i*ceiling(nbMarkers/nbPlot))]) +
                      RotatedAxis() +
                      CenterTitle() +
                      ggtitle(paste0(title, " - ", i, "/", (nbMarkers%/%40)+1, "\n", SO@project.name)) +
                      scale_colour_gradientn(colours = myPalette(100)))
        }
    }
}


createDimPlot <- function(SO, title=NA, cols=NULL){
    ## DimPlot
    p <- DimPlot(SO, reduction = "umap2d", pt.size = 0.8,
                 label = F, cols = cols) +
        ggtitle(title) +
        CenterTitle() +
        NoAxes()
    LabelClusters(p, id = "ident",  fontface = "bold", size = 5)
}

    
doFeaturePlot <- function(SO, geneList, split.by=NULL){
    nbMarkers <- length(geneList)
    geneList <- sort(geneList)
    
    if(is.null(split.by)){
        if (nbMarkers %/% 12 == 0) {
            FeaturePlot(SO, features = geneList,
                        reduction = "umap2d", min.cutoff = "q1", max.cutoff = "q99",
                        order = TRUE, split.by = split.by, label = T)
        } else {
            nbPlot <- (nbMarkers%/%12)+1
            for (i in seq(1, nbPlot)){
                print(FeaturePlot(SO, features = geneList[seq(1+(i-1)*(ceiling(nbMarkers/nbPlot)), i*ceiling(nbMarkers/nbPlot))],
                                  reduction = "umap2d", min.cutoff = "q1", max.cutoff = "q99",
                                  order = TRUE, split.by = split.by, label = T))
            }
        }
        
    } else {
        for (gene in geneList){
            plot1 <- FeaturePlot(SO, features = gene,
                                 reduction = "umap2d", min.cutoff = "q1", max.cutoff = "q99",
                                 order = TRUE, label = T)
            plot2 <- FeaturePlot(SO, features = gene,
                                 reduction = "umap2d", min.cutoff = "q1", max.cutoff = "q99",
                                 order = TRUE, split.by = split.by, label = T) +
                theme(title = element_blank())
            
            plot3 <- ggplot() + theme_void() # empty plot
            
            print((plot3 | plot1| plot3)/plot2)
        }
    }
}

    
createRidgePlot <- function(SO, genesList){
    ## RidgePlot
    for (i in seq(1, length(genesList), 4)){
        print(RidgePlot(SO, features = genesList[seq(i, i+3)], cols=alpha(colors.stage,0.66), ncol = 2))
    }
}


createHeatmapAnn <- function(SO, day = NULL){
    if (is.null(day)){
        title <- "Counts of atlas transfert identity labels\nfrom single datasets to combined dataset"
    } else {
        title <- paste0("Counts of atlas transfert identity labels from ", day, " dataset\nto the sames cells in the combined dataset")
    }
    # create heatmap
    Idents(SO) <- SO@meta.data$day
    transition <- table(SO@meta.data[WhichCells(SO, idents = day), 'orig.celltype'], SO@meta.data[WhichCells(SO, idents = day), 'celltype'])
    dt <- as.data.frame(transition, row.names = names(SO@meta.data[WhichCells(SO, idents = day), 'celltype']))
    dt[dt == 0] <- NA
    names(dt) <- c("singleAnnotation", "combineAnnotation", "counts")
    ggplot(dt, aes(singleAnnotation, combineAnnotation, fill=counts)) +
        geom_tile() +
        geom_text(aes(label = counts), color=if (is.null(day)) "white" else "black", size = 2.5, fontface="bold") +
        #scale_color_manual(values = c(NA, "black"), na.value = NA) +
        coord_fixed() +
        scale_fill_gradient(low = "white", high = "blue", na.value="white") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        guides(fill = guide_colourbar(barwidth = 0.5,
                                      barheight = 20,
                                      title = "Counts")) +
        ggtitle(title) +
        CenterTitle() +
        labs(x = "Celltype annotation on the single datasets", y = paste0("Celltype annotation on the ", SO@project.name, " dataset"))
}
