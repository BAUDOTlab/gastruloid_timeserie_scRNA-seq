extract_topn <- function(markers, topn=20){

    # Rearrange the columns to be more intuitive
    # Order values to be more intuitive
    markers <- markers %>%
        dplyr::arrange(cluster, p_val_adj)
    # Extract top markers per cluster
    topnPos <- markers %>%
        group_by(cluster) %>%
        top_n(n = topn,
              wt = avg_log2FC)
    topnNeg <- markers %>%
        group_by(cluster) %>%
        top_n(n = -topn,
              wt = avg_log2FC)
              
    topnMarkers <- bind_rows(topnPos, topnNeg)
    topnMarkers <-  topnMarkers %>%
        dplyr::arrange(cluster, p_val_adj)
    
    return(topnMarkers)
}
