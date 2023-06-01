#!/usr/bin/env Rscript

library(optparse)
library(purrr)
library(callr)

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(purrr))

option_list <- list(
    make_option(c("-N", "--analysis_name"), action="store", default=NA, type='character',
                help="name of the analysis. It might be already existing, knowing that
                it might overwrite existing files. This argument is MANDATORY for multiple
                time points analysis, with the name of analysis in which multiple single
                datasets were analyzed"),
    
    make_option(c("--new_analysis"), action="store_true", default=FALSE, type='logical',
                help="whether to start a new bunch of analysis. When FALSE, it will overwrite
                the last analysis if already done, if TRUE, a new output directory will be created
                to store all the daily dataset analysis and combined dataset analysis"),
    
    make_option(c("-i", "--input_dataset"), action="store", default=NA, type='character',
                help=""),
    
    make_option(c("--merge"), action="store_true", default=FALSE, type='logical',
                help="if running multiple time points analysis, whether to merge the data,
                otherwise, data will be integrated with Seurat integration process"),
    
#     make_option(c("-l", "--mito_low"), action="store", default=1.5, type='numeric',
#                 help="lower threshold of mitochondrial expression percentage to filter out", metavar = "[0:99]"),
    
    make_option(c("-m", "--mito_high"), action="store", default=10, type='numeric',
                help="", metavar = "[1:100]"),
    
    make_option(c("-q", "--ribo_low"), action="store", default=25, type='numeric',
                help=""),
    
#     make_option(c("-r", "--ribo_high"), action="store", default=45, type='numeric',
#                 help=""),
    
    make_option(c("-f", "--min_feat"), action="store", default=200, type='integer',
                help="minimum number of features detected in every cells"),
    
    make_option(c("-c", "--min_cells"), action="store", default=3, type='integer',
                help="minimum number in which a feature must be detected"),
    
    make_option(c("--min_counts"), action="store", default=1500, type='integer',
                help=""),
    
    make_option(c("--max_counts"), action="store", default=150000, type='integer',
                help=""),
    
    make_option(c("-n", "--norm_method"), action="store", default="LogNormalize", type='character',
                help="normalization method. Imply to not use the same Seurat function in the pipeline.
                The choice is between:
                \t- NormalizeData, with the parameter normalization.method = \"LogNormalize\"
                \t- SCTransform"),
    
    make_option(c("-v", "--hvg_method"), action="store", default="mvp", type='character',
                help="FindVariableFeatures method. The choice is given between:
                \t- vst, what requires to set the number of features to select --hvg_number
                \t- mvp, where the --hvg_number is useless
                \t- disp, that also requires --hvg_number to be set"),
    
    make_option(c("-w", "--hvg_number"), action="store", default=FALSE, type='integer',
                help="number of highly variable genes to use for the downstream analysis.
                Useful when --hvg_method is set to \"vst\" or \"disp\""),
    
    make_option(c("-d", "--do_scaling"), action="store_true", default=FALSE, type='logical',
                help="whether to compute scaling or not. If not mentionned in the command line,
                it will not be done"),
    
    make_option(c("-p", "--pca_npcs"), action="store", default=50, type='integer',
                help="number of PCs to compute for the Seurat function RunPCA"),
    
    make_option(c("--pca_print"), action="store", default=10, type='integer',
                help="number of features to print from the top of each PC"),
    
    make_option(c("-t", "--top_pcs"), action="store", default=30, type='integer',
                help="number of PCs to select for the downstream analysis"),
    
    make_option(c("-o", "--doublet_rate"), action="store", default=8, type='integer',
                help="fraction of doublets estimated by 10XGenomics,
                according to the technical specifiactions"),
    
    make_option(c("-s", "--selected_resolution"), action="store", default=1, type='double',
                help="value of the resolution for the Seurat function FindClusters"),
    
    make_option(c("-a", "--algo_clustering"), action="store", default=4, type='integer',
                help="select the algorithm to run the Seurat function FindClusters
                1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement;
                3 = SLM algorithm; 4 = Leiden algorithm (recommended)"),
                
	make_option(c("-k", "--markers_number"), action="store", default=20, type='integer',
				help="reduce the list of markers to the top [-k] for each cluster at the
				[--selected_resolution]"),
	
	make_option(c("-C", "--combine_meth"), action="store", default=NA, type='character',
				help="indicate the way to combine the datasets. The value has to be one of
				\"merged\", \"blkH\" to perform Harmony integration in one block, \"blkS\" to
				perform Seurat integration in on block and \"seqS\" to perform sequential
				Seurat integration")
)


parsed <- OptionParser(option_list=option_list, add_help_option=TRUE, formatter = IndentedHelpFormatter)
opt = parse_args(parsed)

###########################
# ---- OPTIONS TESTS ---- #
###########################
# input_dataset
if (is.na(opt$input_dataset)){
    callr::rscript("pipeline_launcher.R", cmdargs = "-h")
    cat("You should enter an input dataset among the values:
    - individual dataset analysis
        lab_day_04  lab_day_05
        lab_day_06  lab_day_11
    - time series analysis
        lab_3_days  lab_4_days
        ")
    quit(status = 1)
} else if (opt$input_dataset %in% c("lab_day_04", "lab_day_05", "lab_day_06", "lab_day_11")){
    # TEST POUR L'INTEGRATION DES DONNEES
    if (opt$merge){
        callr::rscript("pipeline_launcher.R", cmdargs = "-h")
        cat("Data of a single time point cannot be merged
            ")
        quit(status = 1)
    } else {
        opt$do_integ <- FALSE
        opt$do_merge <- FALSE
    }
} else if (opt$input_dataset %in% c("lab_3_days", "lab_4_days")){
    if (opt$new_analysis){
        callr::rscript("pipeline_launcher.R", cmdargs = "-h")
        cat("A combined analysis of multiple time points can't be a new
        analysis. Daily datasets must be analyzed beforehand
        ")
        quit(status = 1)
    }
    cat("You are running a multiple time points analysis. The daily datasets have
    to be analyzed beforehand.
    The next parameters will not be taken into account
    - at all:
        --mito_low, --mito_high, --ribo_low, --ribo_high, --doublet_rate
    - on your data, but only on the atlas data:
        --min_feat, --min_cells, --min_counts, --max_counts
    ")
    if (opt$merge){
        opt$do_integ <- FALSE
        opt$do_merge <- TRUE
        opt$combine_meth <- "merged"
    } else {
        opt$do_integ <- TRUE
        opt$do_merge <- FALSE
    }
} else {
    callr::rscript("pipeline_launcher.R", cmdargs = "-h")
    cat("You should enter an input dataset among the values:
    - individual dataset analysis
        lab_day_04  lab_day_05
        lab_day_06  lab_day_11
    - time series analysis
        lab_3_days  lab_4_days
    ")
    quit(status = 1)
}

# norm_method
if (opt$norm_method == "SCTransform"){
    cat("The SCTransform method is not yet implemented, please use 'LogNormalize'
        ")
    quit(status = 1)
} else if (opt$norm_method != "LogNormalize"){
    callr::rscript("pipeline_launcher.R", cmdargs = "-h")
    cat("The 'LogNormalize' method is the only one implemented as for now, please
    let the default parameter for --norm_method or set it at 'LogNormalize'
        ")
    quit(status = 1)
}

# hvg_method and hvg_number
if (opt$hvg_method != "mvp"){
    if (opt$hvg_method %in% c("vst", "disp")){
        if (opt$hvg_number == 0){
            opt$hvg_number <- 2000
        }
    } else {
        callr::rscript("pipeline_launcher.R", cmdargs = "-h")
        quit(status = 1)
    }
} else if (opt$hvg_number){
    callr::rscript("pipeline_launcher.R", cmdargs = "-h")
    quit(status = 1)
} else if (opt$hvg_number == 0){
    opt$hvg_number <- FALSE
}

# algo_clustering
if (!(opt$algo_clustering %in% seq(1, 4))){
    callr::rscript("pipeline_launcher.R", cmdargs = "-h")
    cat("
    The choice of the clustering algorithm is done by using integers
    from 1 to 4. Please see above documentation --algo_clustering.
        ")
    quit(status = 1)
}


param.list <- list(
    analysis.name = opt$analysis_name,
    new.analysis = opt$new_analysis,
    dataset = opt$input_dataset,
    do.integ = opt$do_integ,
    do.merge = opt$do_merge,
#     mito.low = opt$mito_low,
    mito.high = opt$mito_high,
    ribo.low = opt$ribo_low,
#     ribo.high = opt$ribo_high,
    min.feat = opt$min_feat,
    min.cells = opt$min_cells,
    min.counts = opt$min_counts,
    max.counts = opt$max_counts,
    norm.meth = opt$norm_method,
    hvg.meth = opt$hvg_method,
    hvg.num = opt$hvg_number,
    do.scale = opt$do_scaling,
    pca.npcs = opt$pca_npcs,
    pca.print = opt$pca_print,
    top.pcs = opt$top_pcs,
    dblt.rate = opt$doublet_rate,
    res = opt$selected_resolution,
    algo.cluster = opt$algo_clustering,
    top.markers = opt$markers_number,
    combine.meth = opt$combine_meth
)


#basePath <- "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis" # local server
basePath <- "/shared/projects/mothard_in_silico_modeling/seurat_analysis" # IFB server

# multi TP analysis ?
if (opt$input_dataset %in% c("lab_3_days", "lab_4_days")){
    # YES
    
    # Is it a former analysis ? Did you named it ?
    if (!opt$new_analysis && !is.na(opt$analysis_name)){
        # YES and YES
        
        # Does the named analysis already exist ?
        if (length(grep(opt$analysis_name, list.dirs(basePath, full.names = FALSE, recursive = FALSE), value = TRUE)) == 1){
            # YES
            
            combination <- if (opt$do_integ) "integration" else "merge"
            rmarkdown::render("combinedDatasets_pipeline.Rmd", params = param.list, output_file = paste0("combinedData_", opt$analysis_name, "_merged_", opt$input_dataset, ".html"))
        } else {
            # NO
            
            callr::rscript("pipeline_launcher.R", cmdargs = "-h")
            cat("A multi time point analysis MUST be named (--analysis_name), and that
            name is the one where independently time point data where already preprocessed
        ")
            quit(status = 1)
        }
    } else {
        # any combination YES and NO, NO and YES, NO and NO
        
        callr::rscript("pipeline_launcher.R", cmdargs = "-h")
        cat("A multi time point analysis MUST be named (--analysis_name), and that
            name is the one where independently time point data where already preprocessed
        ")
    }
} else {
    # NO: single TP analysis
    
    # Did you named the analysis ?
    if (!is.na(opt$analysis_name)){
        # YES
        
        # Is it a new analysis that doesn't exist ?
        # OR
        # Is it a former analysis that already exists ?
        if ((opt$new_analysis && length(grep(opt$analysis_name, list.dirs(basePath, full.names = FALSE, recursive = FALSE), value = TRUE)) == 0) ||
            (!opt$new_analysis && length(grep(opt$analysis_name, list.dirs(basePath, full.names = FALSE, recursive = FALSE), value = TRUE)) == 1)){
            # YES and NO, NO and YES
            
            rmarkdown::render("singleDataset_pipeline.Rmd", params = param.list, output_file = paste0("singleDataset_", opt$analysis_name, "_", opt$input_dataset, ".html"))
            # # select pipeline to run according the analysis_name
            # if (opt$analysis_name == "removedDoublets"){
            #     rmarkdown::render("singleDataset_pipeline.Rmd", params = param.list, output_file = paste0("singleDataset_", opt$analysis_name, "_", opt$input_dataset, ".html"))
            # } else if (opt$analysis_name == "withoutDoublet_subset"){
            #     rmarkdown::render("no_df_subset_pipeline.Rmd", params = param.list, output_file = paste0("singleDataset_", opt$analysis_name, "_", opt$input_dataset, ".html"))
            # }
        } else {
            # incompatible options like:
            # new analysis with a pre-existing analysis name OR
            # not new analysis, but that doesn't exist yet
            
            callr::rscript("pipeline_launcher.R", cmdargs = "-h")
            cat("Incompatible options were set:
            if you call the program with the --analysis_name option, the directory
            should already exists, OR
            call also the --new_analysis to make it as a new analysis, named by
            --analysis_name
                ")
        }
    } else {
        # NO
        
        # run the analysis, and manage automatically the output directory whether
        # it is a new analysis or not, based on a default name "analysis_X"
        rmarkdown::render("singleDataset_pipeline.Rmd", params = param.list, output_file = paste0("singleDataset_", opt$analysis_name, "_", opt$input_dataset, ".html"))
    }
}



# if (opt$input_dataset %in% c("lab_day_04", "lab_day_05", "lab_day_06", "lab_day_11")){
#     rmarkdown::render("singleDataset_pipeline.Rmd", params = param.list, output_file = paste0("singleDataset_", opt$input_dataset, ".html"))
# } else if (opt$input_dataset %in% c("lab_3_ days", "lab_4_days")){
#     combination <- if (opt$do_integ) "integrated" else "merged"
#     rmarkdown::render("combinedDatasets_pipeline.Rmd", params = param.list, output_file = paste0("combinedData_", combination, "_", opt$input_dataset, ".html"))
# }

#rmarkdown::render("no_df_subset_pipeline.Rmd", params = param.list, output_file = paste0("no_df_subset_", opt$input_dataset, ".html"))
#rmarkdown::render("no_df_pipeline.Rmd", params = param.list, output_file = paste0("no_df_", opt$input_dataset, ".html"))


# l <- list(a="1", b=1, c=list(a="1", b=1))
# yaml::write_yaml(l, "list.yaml")


