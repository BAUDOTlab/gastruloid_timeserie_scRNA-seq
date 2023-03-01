#!/opt/anaconda3/envs/urd/bin/R

args <- commandArgs(trailingOnly = TRUE)
#args <- "urd_parameters.csv"
parameterList <- read.csv(args[1])
#parameterList <- parameterList[,1:c(dim(parameterList)[2]-2)]

NBR_SET_OF_PARAMS <- dim(parameterList)[1]          # nbr of line in the input file

for(s in seq(1, NBR_SET_OF_PARAMS)){
    # parameters management
    param.list <- list(
        knn = parameterList[s, "dm.knn"],
        sigma = parameterList[s, "dm.sigma"],
        floodNsim = parameterList[s, "floodPT_n.sim"],
        cellsForward = parameterList[s, "pseudoT_optimal.cells.forward"],
        cellsBack = parameterList[s, "pseudoT_max.cells.back"],
        RWmaxSteps = parameterList[s, "randomWfromTips_max.steps"],
        #treeMeth = parameterList[s, "buildTree_divergence.method"],
        treeMeth = paste0("'", parameterList[s, "buildTree_divergence.method"], "'"),
        cellsPerPTbin = parameterList[s, "buildTree_cells.per.PT.bin"],
        binsPerPTwindow = parameterList[s, "buildTree_bins.per.PT.window"],
        treePtreshold = parameterList[s, "buildTree_p.tresh"]
    )
    
    parNames <- names(param.list)
    param.vect <- NULL
    for(i in 1:length(param.list)){
        param.vect <- append(param.vect, paste(parNames[i], param.list[[i]], sep = " = "))
    }
    param.str <- paste(param.vect, collapse = ", ")
    
    # # writing management
    toWrite <- paste0("Rscript -e \"rmarkdown::render('10.2_parameterOptimization.Rmd', params = list(", param.str,
                      "), output_file = '/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/trajectory/ciml5/paramOptim/reports/knn.", param.list$knn,
                      "_sigma.", param.list$sigma,
                      "_floodNsim.", param.list$floodNsim,
                      "_cellsForward.", param.list$cellsForward,
                      "_cellsBack.", param.list$cellsBack,
                      "_RWmaxSteps.", param.list$RWmaxSteps,
                      "_treeMeth.", parameterList[s, "buildTree_divergence.method"],
                      "_cellsPerPTbin.", param.list$cellsPerPTbin,
                      "_binsPerPTwindow.", param.list$binsPerPTwindow,
                      "_treePtreshold.", param.list$treePtreshold,
                      ".html')\""
    )
    
    cat(toWrite, file = "runURD_parallel.txt", sep = "\n", append = TRUE)
}

