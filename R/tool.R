##' @title Perform omics datasets evaluation
##' @description Perform omics datasets evaluation
##' @param data_dir A folder contains multiple omics datasets in tsv format.
##' Row: genes/proteins, column: samples.
##' @param x2 A tsv format file which contains protein or gene expression data.
##' When the value of parameter is not NULL, the data in this file will be used
##' to perform correlation analysis between this dataset and the dataset in
##' data_dir
##' @param sample_list A file in tsv format contains sample class, batch and
##' order information.
##' @param data_type The quantification data type in folder data_dir:
##' protein, gene. Default is protein.
##' @param class_for_cor The class of samples which will be used to perform
##' correlation analysis. A charactar vector.
##' @param class_for_fun The class of samples which will be used for function
##' prediction. A charactar vector.
##' @param class_for_ml The class of samples which will be used for phenotype
##' prediction. A sample list file or a charactar vector.
##' @use_common_features_for_func_pred whether or not to use common
##' protien/genes for function prediction. Default is TRUE.
##' @param class_color The color for class.
##' @param out_dir Output folder
##' @param cpu The number of CPUs. Default is 0 and all available CPUs will be
##' used.
##' @param missing_value_cutoff The cutoff of missing value filtering. Default
##' is 0.5.
##' @param species Default is human.
##' @return The path of HTML report.
##' @author Bo Wen \email{wenbostar@@gmail.com}
run_omics_evaluation=function(data_dir=NULL,x2=NULL,sample_list=NULL,data_type="protein",
                              class_for_cor=NULL,
                              class_for_fun=NULL,
                              class_for_ml=NULL,
                              use_common_features_for_func_pred=TRUE,
                              class_color=NULL,out_dir="./",cpu=0,
                              missing_value_cutoff=0.5,
                              species="human"){
    res <- list()

    ## input parameters
    res$input_parameters <- list()
    res$input_parameters$sample_list_file <- sample_list
    res$input_parameters$sample_list <- read.delim(sample_list,stringsAsFactors = FALSE)
    res$input_parameters$data_dir <- data_dir
    res$input_parameters$data_type <- data_type
    res$input_parameters$class_for_cor <- class_for_cor
    res$input_parameters$class_for_fun <- class_for_fun
    res$input_parameters$class_for_ml <- class_for_ml
    res$input_parameters$missing_value_cutoff <- missing_value_cutoff
    res$input_parameters$x2 <- x2

    ## import data
    input_data_files <- list.files(path = data_dir,pattern = ".tsv",
                                   full.names = TRUE,include.dirs = TRUE)

    ## the number of datasets
    res$input_parameters$n_datasets <- length(input_data_files)
    dataset_names <- basename(input_data_files) %>%
        str_replace_all(pattern = ".tsv",replacement = "")

    x1 <- list()
    for(i in 1:length(input_data_files)){
        x1[[i]] <- import_data(input_data_files[i],sample_list = sample_list,
                               data_type=data_type,
                               missing_value_cutoff = missing_value_cutoff)
        x1[[i]]@ID <- dataset_names[i]
    }
    names(x1) <- dataset_names

    res$input_parameters$datasets <- x1

    ## rna
    dat2 <- NULL
    if(!is.null(x2)){
        dat2 <- import_data(x2,sample_list = sample_list,
                            missing_value_cutoff = missing_value_cutoff)
    }

    save(x1,file = paste(out_dir,"/input_x.rda",sep=""))

    ## run basic metrics
    if(is.null(class_color)){
        class_names <- res$input_parameters$sample_list$class %>% unique
        class_color_data <- data.frame(class = class_names, col=rainbow(length(class_names)))
        class_color <- paste(out_dir,"/class_color.tsv",sep="")
        write_tsv(class_color_data,path=class_color)
    }
    basic_metrics_res <- calc_basic_metrics(x1,class_color=class_color,out_dir=out_dir,cpu=cpu)
    saveRDS(basic_metrics_res,file = paste(out_dir,"/basic_metrics_res.rds",sep=""))
    res$basic_metrics <- basic_metrics_res

    ## batch effect evaluation
    n_batch <- res$input_parameters$sample_list$batch %>% unique() %>% length
    cat("Batch number:",n_batch,"\n")
    if(n_batch >= 2){
        res$batch_effect_metrics <- calc_batch_effect_metrics(x1)
        ##
        res$pca_batch_plot <- plot_pca(basic_metrics_res$datasets,out_dir = out_dir,prefix = "pca_view")
        res$pca_batch_plot_13 <- plot_pca(basic_metrics_res$datasets,out_dir = out_dir,prefix = "pca_view_13",pc="13")
    }

    ## complex
    ## function prediction
    if(species=="human"){
        network_data <- import_network_data(type = data_type)
        network_raw_res <- calc_network_corr(x1,network_data,sample_class=class_for_cor,
                                         missing_value_ratio=0.00001,cpu=cpu)
        network_table_res <- plot_network_cor(network_raw_res$result,out_dir = out_dir,
                                              prefix = "network")
        saveRDS(network_raw_res,file = paste(out_dir,"/network_raw_res.rds",sep=""))
        saveRDS(network_table_res,file = paste(out_dir,"/network_table_res.rds",sep=""))

        res$network_result <- network_raw_res
        res$network_table <- network_table_res
    }

    if((data_type == "protein" || data_type == "gene") && !is.null(x2)){
        ## rna protein
        protein_rna_res <- calc_protein_rna_corr(x1,dat2,sample_class = class_for_cor,
                                                 out_dir = out_dir,cpu = cpu,
                                                 missing_value_ratio=missing_value_cutoff)
        res$protein_rna <- protein_rna_res
        saveRDS(protein_rna_res,file = paste(out_dir,"/protein_rna_res.rda",sep=""))
    }

    ## phenotype prediction
    if(!is.null(class_for_ml)){
        if(!is.null(class_for_ml) && file.exists(class_for_ml)){
            ml_res <- calc_ml_metrics(x1,sample_list=class_for_ml,cpu=cpu)
        }else{
            ml_res <- calc_ml_metrics(x1,sample_class=class_for_ml,cpu=cpu)
        }
        res$ml <- ml_res
    }

    ## function prediction
    if(species=="human"){
        save(x1,missing_value_cutoff,cpu,out_dir,file = "data_for_function_prediction.rda")
        fp_res <- calc_function_prediction_metrics(x1,missing_value_cutoff=missing_value_cutoff,
                                               sample_class=class_for_fun,
                                               use_all=!use_common_features_for_func_pred,
                                               cpu=cpu,out_dir=out_dir,prefix="omicsev")
        res$fun_pred <- fp_res
    }
    ##
    rfile <- paste(out_dir,"/final_res.rds",sep="")
    saveRDS(res,file = rfile)
    report_file <- paste(out_dir %>% normalizePath(),"/final_evaluation_report.html",sep="")
    #run_reporter(rfile,report_file,x2=x2,n_batch = n_batch)
    run_reporter(rfile,report_file)
    return(report_file)
}


import_data=function(file,sample_list=NULL,ratio_pair = NULL,
                     missing_value_cutoff=0.5,
                     data_type="protein"){
    x <- read.delim(file,stringsAsFactors = FALSE,check.names = FALSE)
    dat <- x %>% rename(name=ID)
    para <- new("metaXpara")
    para@rawPeaks <- dat
    para@sampleList <- read.delim(sample_list,stringsAsFactors = FALSE,check.names = FALSE)
    para <- reSetPeaksData(para)
    if(is.null(ratio_pair)){
        para@ratioPairs <- ""
    }else{
        para@ratioPairs <- ratio_pair
    }

    para@dataType <- data_type
    cat("Data type:",para@dataType,"\n")


    #para@missValueImputeMethod <- "none"
    para@missValueImputeMethod <- "knn"
    n_features <- para@peaksData$ID %>% unique() %>% length
    para_m <- metaX::filterPeaks(para,ratio = missing_value_cutoff)
    n_features_m <- para_m@peaksData$ID %>% unique() %>% length
    cat("file:",file,", total features:",n_features," 50% missing:",n_features_m,"\n")
    return(para)
}

get_metrics=function(x, metric = "total_features"){
    res <- sapply(x,function(y)y[[metric]])
    return(res)
}


run_reporter=function(x,out_file="test.html",x2=NULL,n_batch=0){
    # if(is.null(x2)){
    #     if(n_batch >= 2){
    #         rmd <- system.file("report/report_noX2.rmd",package = "OmicsEV")
    #     }else{
    #         rmd <- system.file("report/report_noX2_nobatch.rmd",package = "OmicsEV")
    #     }
    # }else{
    #     if(n_batch >= 2){
    #         rmd <- system.file("report/report.rmd",package = "OmicsEV")
    #     }else{
    #         rmd <- system.file("report/report_nobatch.rmd",package = "OmicsEV")
    #     }
    # }
    rmd <- system.file("report/report.rmd",package = "OmicsEV")
    x <- normalizePath(x)
    work_dir <- dirname(x)
    cat("Input file:",x,"\n")
    cat("Work dir:",work_dir,"\n")
    render(rmd,params = list(input=x),output_file = out_file, knit_root_dir=work_dir)
}


launch_gui=function(){
    runApp(appDir = system.file("app",package = "OmicsEV"))
}





