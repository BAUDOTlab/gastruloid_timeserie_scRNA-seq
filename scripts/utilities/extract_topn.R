extract_topn <- function(markers, topn=20){

    # Rearrange the columns to be more intuitive
    # Order values to be more intuitive
    markers <- markers %>%
        dplyr::arrange(cluster, p_val_adj)
    # Extract top markers per cluster
    topnMarkers <- markers %>%
        group_by(cluster) %>%
        top_n(n = topn,
              wt = avg_log2FC)
    
    return(topnMarkers)
}
