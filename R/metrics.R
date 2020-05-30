
#
import_network_data=function(file=NULL,type="protein",species="human"){
    if(species=="human"){
        if((type=="gene" || type=="protein") && is.null(file)){
            x <- readRDS(system.file("extdata/corum.rds",package = "OmicsEV"))
        }else if(type=="phosphopeptide"){
            x <- readRDS(system.file("extdata/psp.rds",package = "OmicsEV"))
        }
    }else if(species=="drosophila"){
        if((type=="gene" || type=="protein") && is.null(file)){
            x <- readRDS(system.file("extdata/flybase_ppi.rds",package = "OmicsEV"))
        }
    }else{
        cat("species:",species,"\n")
        stop("The species is not supported!")
        x <- NULL
    }
    return(x)
}

import_network_data_from_table=function(file){
    a <- readr::read_tsv(file)
    x <- a$complex
    n_p <- sapply(x,function(y){
        a=unique(unlist(strsplit(split=";",x=y)));
        length(a)})
    # filter complex only having one protein
    x <- x[n_p>=2]
    x <- sapply(x, function(y){
        a=unique(unlist(strsplit(split=";",x=y)));
        a=sort(a);paste(a,sep=";",collapse = ";")}) %>%
        unique
    aa <- sapply(x, function(y){
        t(combn(unlist(strsplit(split = ";",x=y)),2))})
    names(aa) <- NULL
    aa <- do.call(rbind,aa)
    res <- list()
    res$complex <- apply(aa,1,function(xx){
        sort(str_replace_all(xx,pattern=" ",replacement=""))}) %>%
        t %>%
        dplyr::as_data_frame() %>%
        distinct() %>%
        filter(V1!="",V2!="")
    res$complex <- res$complex %>% filter(V1!=V2)

    ## random complex
    rnd <- c(res$complex$V1,res$complex$V2) %>%
        unique %>%
        sort %>%
        combn(m=2) %>%
        t %>%
        dplyr::as_data_frame()

    yyy <- paste(res$complex$V1,res$complex$V2,sep="|")
    rnd_p <- paste(rnd$V1,rnd$V2,sep="|")
    rnd <- rnd[ !(rnd_p %in% yyy), ]
    res$rnd_complex <- rnd
    return(res)
}


import_data_from_corum=function(file){
    a <- read.csv(file,stringsAsFactors = FALSE)
    x <- a$subunits.Gene.name.
    n_p <- sapply(x,function(y){
        a=unique(unlist(strsplit(split=";",x=y)));
        length(a)})
    # filter complex only having one protein
    x <- x[n_p>=2]
    x <- sapply(x, function(y){
        a=unique(unlist(strsplit(split=";",x=y)));
        a=sort(a);paste(a,sep=";",collapse = ";")}) %>%
        unique
    aa <- sapply(x, function(y){
        t(combn(unlist(strsplit(split = ";",x=y)),2))})
    names(aa) <- NULL
    aa <- do.call(rbind,aa)
    res <- list()
    res$complex <- apply(aa,1,function(xx){
        sort(str_replace_all(xx,pattern=" ",replacement=""))}) %>%
        t %>%
        dplyr::as_data_frame() %>%
        distinct() %>%
        filter(V1!="",V2!="")
    res$complex <- res$complex %>% filter(V1!=V2)

    ## random complex
    rnd <- c(res$complex$V1,res$complex$V2) %>%
        unique %>%
        sort %>%
        combn(m=2) %>%
        t %>%
        dplyr::as_data_frame()

    yyy <- paste(res$complex$V1,res$complex$V2,sep="|")
    rnd_p <- paste(rnd$V1,rnd$V2,sep="|")
    rnd <- rnd[ !(rnd_p %in% yyy), ]
    res$rnd_complex <- rnd
    return(res)
}

# wget http://www.droidb.org/data/DroID_v2018_08/flybase_ppi.txt
# file = "inst/extdata/flybase_ppi.txt"
# import_data_from_flybase_ppi(file)
import_data_from_flybase_ppi=function(file){

    a <- read.delim(file,stringsAsFactors = FALSE)
    x <- a %>% select(SYMBOL1,SYMBOL2)
    names(x) <- c("V1","V2")
    x <- x %>% filter(V1!=V2) %>% distinct() %>%
        apply(1,function(y) {sort(y)}) %>%
        t %>%
        as.data.frame(stringsAsFactors=FALSE) %>%
        distinct()
    all_genes <- sort(unique(c(x$V1,x$V2)))

    res <- list()
    res$network <- x

    ## random complex
    rnd <- sort(all_genes) %>%
        combn(m=2) %>%
        t %>%
        as.data.frame(stringsAsFactors=FALSE)

    yyy <- paste(res$network$V1,res$network$V2,sep="|")
    rnd_p <- paste(rnd$V1,rnd$V2,sep="|")
    rnd <- rnd[ !(rnd_p %in% yyy), ]
    res$rnd_network <- rnd
    return(res)
}


calc_network_corr=function(x,y,sample_class=NULL,missing_value_ratio=0.00001,cpu=0){

    res <- list()
    res$use_sample_class <- sample_class

    if(!is.null(sample_class)){
        cat("Only use samples from class: ",
            paste(sample_class,collapse = ","),"\n")
        x <- lapply(x, function(xx){
                rm_samples <- xx@sampleList %>% filter(!(class %in% sample_class))
                yy <- removeSample(xx,rsamples = rm_samples$sample)
                cat("Use samples:\n")
                out <- yy@peaksData %>% dplyr::select(sample,class) %>%
                    dplyr::distinct() %>%
                    dplyr::group_by(class) %>%
                    dplyr::summarise(n=n())
                print(out)
                return(yy)
        })
    }else{
        cat("Use all samples!\n")
    }

    res$missing_value_ratio <- missing_value_ratio

    if(missing_value_ratio>0){
        cat("Filter missing value with cutoff: ",missing_value_ratio,"\n")
        x <- lapply(x, function(xx){
            yy <- xx %>% filterPeaks(ratio=missing_value_ratio)
            return(yy)
        })
    }else{
        cat("No missing value filtering!\n")
    }


    for(i in 1:length(x)){
        if(i == 1){
            common_ids <- x[[1]]@peaksData$ID
        }else{
            common_ids <- intersect(common_ids,x[[i]]@peaksData$ID)
        }
    }


    common_ids <- unique(common_ids)
    res$common_ids <- common_ids
    res$n_used_pairs <- nrow(y$network)

    cat("Common genes:",length(common_ids),"\n")
    cat("Used complexes:",nrow(y$network),"\n")
    if(nrow(y$rnd_network) > 10*nrow(y$network)){
        set.seed(1900)
        n_rnd <- nrow(y$rnd_network)
        n <- nrow(y$network)
        ind <- sample(1:n_rnd,2*n)
        y$rnd_network <- y$rnd_network[ind,]

        cat("Total random complexes:",n_rnd,"\n")
        cat("Used random complexes:",nrow(y$rnd_network),"\n")
    }

    if(cpu==0){
        cpu <- detectCores()
    }
    if(cpu > length(x)){
        cpu <- length(x)
    }
    cat("Used cpus:",cpu,"\n")
    cl <- makeCluster(getOption("cl.cores", cpu))
    clusterExport(cl, c("calc_pair_cor"),envir=environment())
    clusterExport(cl, c("calcCor"),envir=environment())
    #clusterEvalQ(cl,library("dplyr"))

    tmp_res <- parLapply(cl,x,fun = calc_pair_cor,pdat=y,idset=common_ids)
    res$result <- tmp_res

    stopCluster(cl)

    return(res)
}


calc_pair_cor=function(para,pdat,idset=NULL){
    ## filter genes
    valid_ids <- para@peaksData$ID %>% unique
    pdat$network <- pdat$network %>%
        filter(V1 %in% valid_ids, V2 %in% valid_ids)
    pdat$rnd_network <- pdat$rnd_network %>%
        filter(V1 %in% valid_ids, V2 %in% valid_ids)

    ## only consider the genes provided by users
    if(!is.null(idset)){
        valid_ids <- idset
        pdat$network <- pdat$network %>%
            filter(V1 %in% valid_ids, V2 %in% valid_ids)
        pdat$rnd_network <- pdat$rnd_network %>%
            filter(V1 %in% valid_ids, V2 %in% valid_ids)
    }

    pdat$network <- pdat$network %>% mutate(nid=1:nrow(pdat$network))
    dat <- pdat$network %>% gather(key="network",value="ID",-nid)
    m <- inner_join(dat,para@peaksData,by="ID")
    res <- list()
    res$network <- m %>% group_by(nid) %>%
        do(cor=calcCor(.)) %>%
        ungroup() %>%
        as.data.frame() %>%
        mutate(cor=unlist(cor))

    ## random complex
    pdat$rnd_network <- pdat$rnd_network %>% mutate(nid=1:nrow(pdat$rnd_network))
    dat <- pdat$rnd_network %>% gather(key="network",value="ID",-nid)
    m <- inner_join(dat,para@peaksData,by="ID")
    res$rnd_network <- m %>% group_by(nid) %>%
        do(cor=calcCor(.)) %>%
        ungroup() %>%
        as.data.frame() %>%
        mutate(cor=unlist(cor))
    return(res)
}

calcCor=function(x){
    x <- x %>% dplyr::select(ID,sample,value,nid) %>%
        spread(key=ID,value=value) %>%
        as.data.frame()
    res <- cor(x[,3],x[,4],use = "com",method = "spearman")
    return(res)
}


plot_network_cor=function(x, out_dir="./",prefix="test"){

    res <- list()

    dataset_name <- names(x)
    for(i in 1:length(dataset_name)){
        x[[i]]$network$nid <- dataset_name[i]
        x[[i]]$rnd_network$nid <- dataset_name[i]
    }

    ks_res <- x %>% lapply(function(x){
        x1 <- unlist(x$network$cor)
        x2 <- unlist(x$rnd_network$cor)
        rr <- ks.test(x1,x2)
        rr$statistic
    })
    res$ks <- ks_res
    #ks_file <- paste(out_dir,"/",prefix,"-ks.txt",sep="")
    #write(ks_res,file = ks_file)
    #ppi_cor_res$res11$complex$nid <- "RNA"

    network_data_list <- lapply(x, function(y){
        return(y$network)
    })

    network_data_df <- bind_rows(network_data_list) %>% mutate(type="IntraComplex")

    rnd_network_data_list <- lapply(x, function(y){
        return(y$rnd_network)
    })

    rnd_network_data_df <- bind_rows(rnd_network_data_list) %>% mutate(type="InterComplex")


    dat <- rbind(network_data_df,rnd_network_data_df)
    save(dat,file = "dat.rda")
    fig <- paste(out_dir,"/",prefix,"-complex_all_boxplot.png",sep="")
    res$network_boxplot <- fig
    png(fig,width = 1000,height = 500,res=160)
    aaaa <- dat %>% group_by(nid) %>%
        filter(type=="InterComplex") %>%
        summarize(c=median(cor))
        #arrange(c)
    dat -> xxxx
    dat$nid <- factor(dat$nid,levels = aaaa$nid )
    gg <- ggplot(dat,aes(x=nid,y=cor,colour=type))+
        #geom_violin(position="dodge")+
        geom_boxplot(position="dodge")+
        xlab("Dataset")+
        ylab("Correlation")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(gg)
    dev.off()
    save(dat,file="complex_res.rda")

    cor_res <- dat %>% group_by(nid,type) %>%
        summarise(cor=median(cor)) %>%
        ungroup() %>%
        spread(key="type",value="cor") %>%
        mutate(diff=IntraComplex-InterComplex) %>%
        rename(dataSet=nid)

    ks_df <- data.frame(dataSet=names(ks_res),ks=ks_res %>% unlist())
    res$cor <- merge(cor_res,ks_df,sort=FALSE)


    ## save plot data
    fig_data_file <- paste(out_dir,"/",prefix,"-complex_all_boxplot.rds",sep="")
    fig_data <- list()
    fig_data$plot_object <- gg
    fig_data$plot_data <- dat
    saveRDS(fig_data,file = fig_data_file)


    return(res)
}


## protein rna correlation
# "pearson", "kendall", "spearman"
# This function is used to calculate the correlation between protein and mRNA.
calcCorBetweenProteinAndRNA=function(para1,para2,log=TRUE,select_by=1,top_n=1000,
                                     geneset = NULL,
                                     outdir = "./",prefix = "test",use_class=NULL,
                                     reg=FALSE,
                                     cor_method="spearman",valueID="value"){

    if(!is.null(use_class)){
        cat("Only use class ",use_class,"\n")
        if(reg == FALSE){
            para1@peaksData <- para1@peaksData %>% dplyr::filter(class==use_class)
            para2@peaksData <- para2@peaksData %>% dplyr::filter(class==use_class)
        }else{
            para1@peaksData <- para1@peaksData %>% dplyr::filter(grepl(pattern = use_class,x = class))
            para2@peaksData <- para2@peaksData %>% dplyr::filter(grepl(pattern = use_class,x = class))
        }
        cat("Used samples:",length(unique(para1@peaksData$sample)),length(unique(para2@peaksData$sample)),"\n")
    }else{
        cat("Use all samples!\n")
    }

    if(!is.null(geneset)){
        cat("Only use specified gene set\n")
        para1@peaksData <- para1@peaksData %>% dplyr::filter(ID %in% geneset)
        para2@peaksData <- para2@peaksData %>% dplyr::filter(ID %in% geneset)
    }

    x1 <- para1@peaksData %>% dplyr::select(ID,sample,!!valueID)
    x2 <- para2@peaksData %>% dplyr::select(ID,sample,!!valueID)
    m <- merge(x1,x2,by=c("ID","sample"))

    ## firstly, plot a scatterplot: x axis is sd, y axis is cor
    names(m)[3:4] <- c("x","y")
    m$x[m$x <= 0] <- NA
    m$y[m$y <= 0] <- NA
    if(log==TRUE){
        m$x <- log2(m$x)
        m$y <- log2(m$y)
    }

    ## filter genes with all values are missing
    na_row <- apply(m[,3:4],1,function(x){any(is.na(x))})
    m <- m[!na_row,]

    ## Sample wise correlation
    save(m,file = "m.rda")
    sample_wise_cor <- m %>% group_by(sample) %>%
        dplyr::summarise(cor=cor(x,y,use = "com",method = cor_method)) %>%
        dplyr::mutate(label=prefix)

    ## Gene wise correlation
    res <- m %>% group_by(ID) %>% dplyr::summarise(sd_x=sd(x,na.rm = TRUE),
                                                   sd_y=sd(y,na.rm=TRUE),
                                                   cor=cor(x,y,use = "com",method = cor_method))
    #save(para1,para2,valueID,na_row,m,res,file="test111.rda")

    png(paste(outdir,"/",prefix,"-cor-sd.png",sep=""),width = 600,height = 1200,res=120)
    cor_name <- paste("Cor (",cor_method,")",sep="")
    par(mfrow=c(4,2),mar=c(3,3,2,1),mgp=c(1.6,0.6,0))
    plot(res$sd_x,res$cor, xlab="SD (1)",ylab=cor_name,col=rgb(0.5,0.4,0.7,0.4),pch=15,cex=0.6)
    plot(res$sd_y,res$cor, xlab="SD (2)",ylab=cor_name,col=rgb(0.5,0.4,0.7,0.4),pch=15,cex=0.6)
    plot(res$sd_x,res$sd_y, xlab="SD (1)",ylab="SD (2)",col=rgb(0.5,0.4,0.7,0.4),pch=15,cex=0.6)
    hist(res$sd_x,nclass=50,xlab="SD (1)",main="")
    hist(res$sd_y,nclass=50,xlab="SD (2)",main="")
    cor_mean <- sprintf("%.4f",mean(res$cor,na.rm=TRUE))
    cor_md <- sprintf("%.4f",median(res$cor,na.rm=TRUE))
    hist(res$cor,nclass=50,xlab=cor_name,main=paste("Data point: ",nrow(res),"\nmean = ",cor_mean,", median = ",cor_md,sep = ""))
    dev.off()

    rr <- data.frame(n=nrow(res),
                     n5=sum(res$cor>=0.5,na.rm = TRUE),
                     n6=sum(res$cor>=0.6,na.rm = TRUE),
                     n7=sum(res$cor>=0.7,na.rm = TRUE),
                     n8=sum(res$cor>=0.8,na.rm = TRUE),
                     median_cor=median(res$cor,na.rm=TRUE),
                     stringsAsFactors = FALSE)

    fres <- list()
    fres$feature_wise_cor_data <- res
    fres$feature_wise <- rr
    fres$sample_wise <- sample_wise_cor

    return(fres)
}

calc_protein_rna_corr=function(x,rna,sample_class=NULL,out_dir="./",cpu=0,
                               missing_value_ratio=0.5){

    x <- lapply(x, function(xx){
        yy <- metaX::filterPeaks(xx,ratio = missing_value_ratio)
        return(yy)
    })

    if(!is.null(sample_class)){
        cat("Only use samples from class: ",
            paste(sample_class,collapse = ","),"\n")
        x <- lapply(x, function(xx){
            rm_samples <- xx@sampleList %>% filter(!(class %in% sample_class))
            yy <- removeSample(xx,rsamples = rm_samples)
            return(yy)
        })
    }else{
        cat("Use all samples!\n")
    }

    for(i in 1:length(x)){
        if(i == 1){
            common_ids <- x[[1]]@peaksData$ID
        }else{
            common_ids <- intersect(common_ids,x[[i]]@peaksData$ID)
        }
    }

    common_ids <- intersect(common_ids,rna@peaksData$ID)

    common_ids <- unique(common_ids)
    cat("Common IDs:", length(common_ids),"\n")

    if(cpu==0){
        cpu <- detectCores()
    }
    if(cpu > length(x)){
        cpu <- length(x)
    }
    cat("Use cpus:",cpu,"\n")
    cl <- makeCluster(getOption("cl.cores", cpu))

    res <- parLapply(cl,x,fun = calcCorBetweenProteinAndRNA,
                     para2=rna,
                     log=TRUE,
                     select_by=1,
                     top_n=1000,
                     cor_method="spearman",
                     valueID="value")

    stopCluster(cl)

    dataset_name <- names(res)
    for(i in 1:length(dataset_name)){
        res[[i]]$feature_wise$dataSet <- dataset_name[i]
        res[[i]]$sample_wise$dataSet <- dataset_name[i]
        res[[i]]$feature_wise_cor_data$dataSet <- dataset_name[i]
    }

    ## feature wise
    cor_res <- lapply(res, function(y){y$feature_wise}) %>% bind_rows()
    cor_res <- bind_cols(cor_res %>% select(dataSet),cor_res %>% select(-dataSet))


    ## sample wise
    sample_wise_cor_res <- lapply(res, function(y){y$sample_wise}) %>% bind_rows()
    sample_wise_cor_res_table <- sample_wise_cor_res %>% group_by(dataSet) %>%
        summarise(median_cor=median(cor,na.rm = TRUE))

    fres <- list(feature_wise_cor_table = cor_res,
                 sample_wise_cor_table = sample_wise_cor_res_table,
                 data=res,
                 use_sample_class = sample_class)

    ## plot
    gg <- ggplot(sample_wise_cor_res,aes(x=dataSet,y=cor,color=dataSet)) +
        geom_violin()+
        geom_boxplot(width=0.1)+
        geom_text(data = sample_wise_cor_res_table, aes(x = dataSet, y = median_cor, label = sprintf("%.4f",median_cor)),
                  size = 5, vjust = -0.2, color = "black",)+
        xlab("Dataset")+
        ylab("Sample wise correlation")+
        theme(legend.position="none")

    sample_wise_cor_fig <- paste(out_dir,"/sample_wise_cor.png",sep="")
    png(sample_wise_cor_fig,width = 800,height = 400,res=120)
    print(gg)
    dev.off()
    fres$sample_wise_cor_fig <- sample_wise_cor_fig

    # generate CDF and boxplot plots for gene-wise correlation analysis
    ## boxplot
    gene_wise_cor_res_all <- lapply(res, function(y){y$feature_wise_cor_data}) %>% bind_rows()
    fres$gene_wise_cor_boxplot_fig <- plot_feature_wise_cor_boxplot(gene_wise_cor_res_all,out_dir = out_dir,prefix = "OmicsEV")
    ## CDF
    fres$gene_wise_cor_cdf_fig <- plot_feature_wise_cor_cdf(gene_wise_cor_res_all,out_dir = out_dir,prefix = "OmicsEV")

    return(fres)

}

plot_feature_wise_cor_cdf=function(x, out_dir="./",prefix="test"){

    fig <- paste(out_dir,"/",prefix,"-feature_wise_cor_cdf.png",sep="")
    png(fig,width = 400,height = 400,res=150)
    gg <- ggplot(x,aes(cor,group=dataSet,color=dataSet))+
        stat_ecdf(geom = "step",size=0.3)+
        theme(legend.position = c(0, 1),
              legend.justification = c(0, 1),
              legend.key = element_blank(),
              legend.background=element_blank())+
        ylab("CDF")+
        xlab("Pearson correlation")
    print(gg)
    dev.off()

    ## save plot data
    fig_data_file <- paste(out_dir,"/",prefix,"-feature_wise_cor_cdf.rds",sep="")
    fig_data <- list()
    fig_data$plot_object <- gg
    fig_data$plot_data <- x
    saveRDS(fig_data,file = fig_data_file)

    return(fig)
}

plot_feature_wise_cor_boxplot=function(x, out_dir="./",prefix="test"){

    median_cor <- x %>% group_by(dataSet) %>% summarise(median_cor=median(cor,na.rm = TRUE)) %>% arrange(desc(median_cor))
    x$dataSet <- factor(x$dataSet,levels = median_cor$dataSet)

    fig <- paste(out_dir,"/",prefix,"-feature_wise_cor_boxplot.png",sep="")
    png(fig,width = 800,height = 400,res=120)
    gg <- ggplot(x,aes(x=dataSet,y=cor))+
        geom_boxplot(width=0.5,outlier.size=0.2)+
        ylab("Pearson correlation")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(gg)
    dev.off()
    return(fig)

}


calc_basic_metrics=function(x,class_color=NULL,out_dir="./",cpu=0){

    if(!is.null(class_color)){
        class_col <- read_tsv(class_color) %>% filter(!is.na(col))
    }else{
        class_col <- NULL
    }

    dataset_name <- names(x)
    xx <- list()
    for(i in 1:length(dataset_name)){
        xx[[i]] <- list(name=dataset_name[i],data=x[[i]])
    }
    names(xx) <- dataset_name

    ## output
    fres <- list()

    plist <- list()
    plist$scale <- "pareto"
    plist$t <- 3
    plist$missValueRatioQC <- 0.5
    plist$missValueRatioSample <- 0.5
    plist$classCol <- class_col
    plist$pcaLabel <- NA
    plist$center <- TRUE

    plist$pngWidth <- 6.5
    plist$pngHeight <- 6.5
    plist$pdfWidth <- 6
    plist$pdfHeight <- 6
    plist$legendRowBatch <- 4
    plist$legendRowClass <- 4

    if(length(dataset_name) >=2){
        fig <- plot_upset(x,out_dir = out_dir, prefix = "total_identification_upset")
        fres$datasets_id_overlap <- fig
    }
    density_fig <- plot_density(x,out_dir = out_dir, prefix = "sample_wise_density")


    if(cpu==0){
        cpu <- detectCores()
    }
    if(cpu > length(dataset_name)){
        cpu <- length(dataset_name)
    }
    cat("Use cpu:",cpu,"\n")

    if(cpu > 1){
        cl <- makeCluster(getOption("cl.cores", cpu))
        clusterExport(cl, c("metaXpipe"),envir=environment())


        res <- parLapply(cl,xx,fun = run_basic_metrics,plist=plist,out_dir=out_dir)
        stopCluster(cl)
    }else{
        res <- lapply(xx,run_basic_metrics,plist=plist,out_dir=out_dir)
    }


    fres$parameters <- plist
    fres$density_plot <- density_fig
    fres$datasets <- res


    return(fres)
}


plot_upset=function(x,out_dir="./",prefix="test"){

    a <- lapply(x, function(y){
        return(y@peaksData$ID %>% unique())
    })

    fig <- paste(out_dir,"./",prefix,".png",sep = "")
    png(fig, width = 9.5, height = 5.5,res=120, units = "in")

    m  <- ComplexHeatmap::make_comb_mat(a)
    ss = ComplexHeatmap::set_size(m)
    cs = ComplexHeatmap::comb_size(m)
    ht <- ComplexHeatmap::UpSet(m,
               set_order = order(ss),
               comb_order = order(comb_degree(m), -cs),
               top_annotation = HeatmapAnnotation(
                   "Intersection size" = anno_barplot(cs,
                                                        ylim = c(0, max(cs)*1.1),
                                                        border = FALSE,
                                                        gp = gpar(fill = "black"),
                                                        height = unit(4, "cm")
                   ),
                   annotation_name_side = "left",
                   annotation_name_rot = 90))


    ht <- ComplexHeatmap::draw(ht)
    od <- ComplexHeatmap::column_order(ht)
    ComplexHeatmap::decorate_annotation("Intersection size", {
        grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(3, "pt"),
                  default.units = "native", just = "bottom", gp = gpar(fontsize = 9))
    })

    dev.off()
    return(fig)
}




run_basic_metrics=function(x,plist,out_dir="./"){

    res <- list()

    para <- x$data
    dat_name <- x$name
    c_out_dir <- paste(out_dir,"/",dat_name,sep="")
    para@outdir <- c_out_dir
    para@prefix <- dat_name

    if(is.null(para@sampleList) || is.na(para@sampleList) ||
       nrow(para@sampleList) ==0){
        sampleList  <- read.delim(para@sampleListFile,stringsAsFactors = FALSE)
    }else{
        sampleList  <- para@sampleList
    }

    checkSampleList(sampleList)

    if(is.null(para@ratioPairs)){
        stop("Please set the value of ratioPairs!")
    }

    ##
    plsdaPara <- new("plsDAPara")
    plsdaPara@scale <- plist$scale
    plsdaPara@t <- plist$t
    plsdaPara@do <- FALSE
    missValueRatioQC <- plist$missValueRatioQC
    missValueRatioSample <- plist$missValueRatioSample
    classCol <- plist$classCol
    pcaLabel <- plist$pcaLabel
    center <- plist$center

    ## directory structure, data (all figures and other files), report.html
    raw_outdir <- para@outdir
    ## all figures and other files
    para@outdir <- paste(para@outdir,"/data",sep="")
    makeDirectory(para)

    if(is.null(para@sampleList) || is.na(para@sampleList) ||
       nrow(para@sampleList) ==0){
        expDesign  <- read.delim(para@sampleListFile,stringsAsFactors = FALSE)
    }else{
        expDesign  <- para@sampleList
    }

    expDesign$class <- as.character(expDesign$class)
    expDesign$class[is.na(expDesign$class)] <- "QC"

    expClass <- expDesign %>% group_by(class) %>% summarise(n=length(sample))

    expBatch <- expDesign %>% group_by(batch) %>% summarise(n=length(sample))

    message(date(),"\treset the peaksData...")
    para <- reSetPeaksData(para = para)

    ##
    res$total_features <- para@peaksData$ID %>% unique() %>% length()

    ## 2. pre-processing, filter noise peaks
    if(hasQC(para)){
        message(date(),"\tfilter peaks in QC sample:>=",missValueRatioQC)
        pre_para <- para
        para <- filterQCPeaks(para = para,ratio = missValueRatioQC)
        rpeak <- length(unique(pre_para@peaksData$ID)) - length(unique(para@peaksData$ID))

        paste("Remove peaks in QC samples with missing
            value greater than ",100*missValueRatioQC, " percent:",rpeak,".",sep="")
        rm(pre_para)

    }
    message(date(),"\tfilter peaks in non-QC sample:>=",missValueRatioSample)
    pre_para <- para
    para <- filterPeaks(para = para,ratio = missValueRatioSample)
    rpeak <- length(unique(pre_para@peaksData$ID)) - length(unique(para@peaksData$ID))

    ##
    res$valid_features <- para@peaksData$ID %>% unique() %>% length()

    paste("Remove peaks in non-QC samples with missing
            value greater than ",100*missValueRatioSample," percent:",rpeak,".")
    rm(pre_para)

    message("The CV distribution before normalization: ")
    printCV(para,valueID = "value")

    ## 3. Quality control
    message("plot peak number distribution...")

    fig <- plotPeakNumber(para,legendRowBatch=10)
    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")

    ##
    res$features_number_distribution <- fig$fig

    message("plot CV distribution...")
    class_sample <- para@peaksData %>% dplyr::select(class,sample) %>% distinct()
    class_2 <- which(table(class_sample$class)>=2) %>% names()
    cv_para <- para
    cv_para@peaksData <- cv_para@peaksData %>% dplyr::filter(class %in% class_2)
    cv_para@sampleList <- cv_para@sampleList %>% dplyr::filter(class %in% class_2)

    fig <- metaX::plotCV(cv_para)
    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")
    ##
    res$cv_distribution <- fig$fig

    message("plot missing value distribution...")
    fig <- plotMissValue(para,height = 3.7,width = 9.15)
    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")

    ##
    res$missing_value_distribution <- fig$fig

    message("plot peak intensity distribution...")
    fig <- plotIntDistr(para,width = 9.15)
    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")

    ##
    res$features_quant_distribution <- fig$fig

    message("plot TIC distribution...")
    fig <- plotPeakSumDist(para,valueID = "value",height = 3.7,width = 5.5)
    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")

    ticTable <- para@peaksData %>% group_by(sample,batch,class) %>%
        dplyr::summarise(tic=sum(value,na.rm = TRUE)) %>%
        group_by(batch) %>%
        dplyr::summarise(n00=quantile(tic)[1],
                         n25=quantile(tic)[2],
                         n50=quantile(tic)[3],
                         n75=quantile(tic)[4],
                         n100=quantile(tic)[5])
    print(ticTable)

    message("plot average intensity distribution...")
    fig <- plotPeakMeanDist(para,valueID = "value")
    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")

    meanTable <- para@peaksData %>% group_by(sample,batch,class) %>%
        dplyr::summarise(meanIntensity=mean(value,na.rm = TRUE)) %>%
        group_by(batch) %>%
        dplyr::summarise(n00=quantile(meanIntensity)[1],
                         n25=quantile(meanIntensity)[2],
                         n50=quantile(meanIntensity)[3],
                         n75=quantile(meanIntensity)[4],
                         n100=quantile(meanIntensity)[5])
    print(meanTable)

    message(date(),"\tmissing value inputation...")
    para <- missingValueImpute(para,cpu=1)
    paste("The missing value were imputed by ",para@missValueImputeMethod,".",sep="")


    save(para,file = "para.rds")
    fig <- plotCorHeatmap(para = para,valueID = "value",anno = TRUE,
                          samples = NULL,
                          height = 6,width = 6,sortBy="batch",cluster = FALSE,
                          classCol=classCol)

    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")
    ##
    res$sample_wise_cor_heatmap <- fig$fig

    if(hasQC(para)){
        message("plot correlation heatmap...")
        #saveRDS(para,file = "para.rda")
        prefix_w <- para@prefix
        para@prefix <- paste(prefix_w,"-qc",sep = "")
        fig <- plotCorHeatmap(para = para,valueID = "value",anno = TRUE,
                              samples = NA, ## QC sample
                              height = 6,width = 6,cluster = FALSE,
                              classCol=classCol)
        #fig1 <- paste("data/",basename(fig$fig),sep="")
        #fig2 <- paste("data/",basename(fig$highfig),sep="")
        para@prefix <- prefix_w
        res$qc_sample_wise_cor_heatmap <- fig$fig

    }

    #save(para,t,file="test.rda")
    save(para,plist,scale,center,pcaLabel,classCol,file=paste(para@outdir,"/","pca.rda",sep=""))
    ppca <- transformation(para,method = plist$t,valueID = "value")
    ppca <- metaX::preProcess(ppca,scale = plist$scale,center = plist$center,
                              valueID = "value")
    fig <- metaX::plotPCA(ppca,valueID = "value",scale = "none",batch = TRUE,
                          rmQC = FALSE,label = pcaLabel,classColor = classCol,
                          legendRowClass = plist$legendRowClass,
                          legendRowBatch = plist$legendRowBatch,
                          pngWidth = plist$pngWidth,
                          pngHeight = plist$pngHeight,
                          pdfWidth = plist$pdfWidth,
                          pdfHeight = plist$pdfHeight)
    prefix_bak <- ppca@prefix
    ppca@prefix <- paste(ppca@prefix,"-nobatch",sep="")
    fig_nobatch <- metaX::plotPCA(ppca,valueID = "value",scale = "none",batch = FALSE,
                                  rmQC = FALSE,label = pcaLabel,classColor = classCol,
                                  legendRowClass = plist$legendRowClass,
                                  legendRowBatch = plist$legendRowBatch,
                                  pngWidth = plist$pngWidth,
                                  pngHeight = plist$pngHeight,
                                  pdfWidth = plist$pdfWidth,
                                  pdfHeight = plist$pdfHeight)

    ##
    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")
    res$pca_with_batch <- fig$fig
    res$pca_plot_data <- fig$plotdata


    ## plot heatmap
    save(ppca,classCol,file = "heatmap.rda")
    fig <- plotHeatMap(ppca,valueID="value",log=FALSE,rmQC=FALSE,
                        #scale="row",
                        #clustering_distance_rows="euclidean",
                        #clustering_distance_cols="euclidean",
                        #clustering_method="ward.D2",
                        seriation = FALSE,
                        classCol=classCol,
                        show_colnames=FALSE)

    ppca@prefix <- prefix_bak


    #fig1 <- paste("data/",basename(fig$fig),sep="")
    #fig2 <- paste("data/",basename(fig$highfig),sep="")
    res$cluster_heatmap <- fig$fig


    ## remove QC
    ppca@prefix <- paste(ppca@prefix,"-noqc",sep="")
    fig <- metaX::plotPCA(ppca,valueID = "value",scale = "none",batch = TRUE,
                          rmQC = TRUE,label = pcaLabel,classColor = classCol,
                          legendRowClass = plist$legendRowClass,
                          legendRowBatch = plist$legendRowBatch,
                          pngWidth = plist$pngWidth,
                          pngHeight = plist$pngHeight,
                          pdfWidth = plist$pdfWidth,
                          pdfHeight = plist$pdfHeight)
    ## no QC, no batch
    ppca@prefix <- paste(ppca@prefix,"-nobatch",sep="")
    fig <- metaX::plotPCA(ppca,valueID = "value",scale = "none",batch = FALSE,
                          rmQC = TRUE,label = pcaLabel,classColor = classCol,
                          legendRowClass = plist$legendRowClass,
                          legendRowBatch = plist$legendRowBatch,
                          pngWidth = plist$pngWidth,
                          pngHeight = plist$pngHeight,
                          pdfWidth = plist$pdfWidth,
                          pdfHeight = plist$pdfHeight)
    #plotPLSDA(para)

    if(!is.null(para@ratioPairs) && para@ratioPairs != ""){
        ## only perform this analysis when a user provides comparison group
        ## information.
        para@peaksData$valueNorm <- para@peaksData$value
        #para <- peakStat(para = para,plsdaPara = plsdaPara,doROC=FALSE,
        #                 pcaLabel=pcaLabel,classColor = pcaColor)

        qf <- paste("data/",para@prefix,"-quant.txt",sep="")
    }

    saveRDS(para,file = paste(para@outdir,"/",para@prefix,"-result.rds",
                                  sep=""))

    message("Print information about the current R session:")
    writeLines(capture.output(sessionInfo()),
               paste(para@outdir,"/",para@prefix,"-sessionInfo.txt",sep=""))

    return(res)
}


get_sample_data=function(x){
    a <- lapply(x, function(y){
        yy <- y@peaksData %>% dplyr::select(sample,class,batch) %>%
            dplyr::distinct() %>%
            dplyr::mutate(dataSet=y@ID)
        yy$sample <- as.character(yy$sample)
        yy$class <- as.character(yy$class)
        return(yy)
    })
    xx <- dplyr::bind_rows(a)
    class_samples <- xx %>% group_by(class,dataSet) %>%
        dplyr::summarise(samples=n()) %>%
        ungroup() %>%
        tidyr::spread(key=dataSet,value=samples)

    return(class_samples)
}

get_identification_summary_table=function(x,format=TRUE){
    x1 <- get_metrics(x$basic_metrics$datasets,metric="total_features")
    x2 <- get_metrics(x$basic_metrics$datasets,metric="valid_features")
    a <- data.frame(dataSet=names(x1),x1=x1,x2=x2,stringsAsFactors = FALSE)
    if(format==TRUE){
        a <- a %>%
            dplyr::mutate(x1 = formattable::color_bar("lightgreen")(x1),
                          x2 = formattable::color_bar("lightgreen")(x2)) #%>%
            #dplyr::mutate(x11 = ifelse(x1 >= max(x1),
            #                           cell_spec(x1, color = "red", bold = T),
            #                           cell_spec(x1, color = "green", italic = T)),
            #              x22 = ifelse(x2 >= max(x2),
            #                           cell_spec(x2, color = "red", bold = T),
            #                           cell_spec(x2, color = "green", italic = T))) %>%
            #rename(x1=NULL,x2=NULL)


    }
    names(a) <- c("dataSet","# proteins (genes)","# proteins (genes) [50%]")
    return(a)
}

get_full_path=function(x){
    for(i in 1:length(x)){
        x[i] <- normalizePath(x[i])
    }
    return(x)
}


run_kbet=function(x,out_dir="./",prefix="test"){
    res <- lapply(x, function(y){
        xx <- metaX:::getPeaksTable(y,value="value")
        png(paste(y@ID,".png",sep=""))

        if(min(table(xx$batch)) <= 10){
            dat <- kBET(df=xx[,-c(1:4)],batch = xx$batch,plot = TRUE,
                        heuristic = FALSE,k0 = 3)
        }else{
            dat <- kBET(df=xx[,-c(1:4)],batch = xx$batch,plot = TRUE)
        }
        ff <- list(kbet=dat,name=y@ID)
        dev.off()
        return(ff)
    })

    res_table <- lapply(res, function(a){
        y <- a$kbet$summary[1,] %>% mutate(dataSet=a$name) %>% select(dataSet,everything())
    })
    res_table <- bind_rows(res_table)

    dd <- lapply(res, function(a){
        y <- a$kbet
        dat <- data.frame(value=c(y$stats$kBET.expected,y$stats$kBET.observed),
                          Test=c(rep("Expected(random)",length(y$stats$kBET.expected)),
                                 rep("Observed(kBET)",length(y$stats$kBET.observed))),
                          dataSet=a$name,
                          stringsAsFactors = FALSE)
    })

    gdat <- bind_rows(dd)
    fig <- paste(out_dir,"/",prefix,"-kBET-boxplot.png",sep="")
    png(fig,width = 800,height = 500,res=130)
    gg <- ggplot(gdat,aes(x=Test,y=value,color=dataSet))+
        geom_boxplot()+
        ylab("Rejection rate")
        #facet_grid(.~dataSet)
    print(gg)
    dev.off()

    fres <- list(kBET_boxplot=fig,table=res_table,kbet=res)
    return(fres)
}


calc_silhouette_width=function(x,transform_method="3",scale_method="pareto",
                               center=TRUE,missing_value_ratio=0.5){
    res <- lapply(x, function(y){
        y <- metaX::filterPeaks(y,ratio = missing_value_ratio)
        y@missValueImputeMethod <- "knn"
        y <- metaX::missingValueImpute(y,valueID = "value",method="knn")
        ppca <- metaX::transformation(y,method = transform_method,valueID = "value")
        ppca <- metaX::preProcess(ppca,scale = scale_method,center = center,
                                  valueID = "value")
        xx <- metaX:::getPeaksTable(y,value="value")
        #data: a matrix (rows: samples, columns: features (genes))
        #batch: vector or factor with batch label of each cell
        pca.data <- prcomp(xx[,-c(1:4)], center=FALSE) #compute PCA representation of the data
        batch.silhouette <- batch_sil(pca.data, xx$batch)
        return(batch.silhouette)
    })
    res <- res %>% unlist()
    return(res)
}


calc_pca_batch_regression=function(x,transform_method=3,scale_method="pareto",
                               center=TRUE,missing_value_ratio=0.5,format=FALSE,
                               top_pc=10){
    res <- lapply(x, function(y){
        y <- metaX::filterPeaks(y,ratio = missing_value_ratio)
        y@missValueImputeMethod <- "knn"
        y <- metaX::missingValueImpute(y,valueID = "value",method="knn")
        ppca <- metaX::transformation(y,method = transform_method,valueID = "value")
        ppca <- metaX::preProcess(ppca,scale = scale_method,center = center,
                                  valueID = "value")
        xx <- metaX:::getPeaksTable(ppca,value="value")
        #data: a matrix (rows: samples, columns: features (genes))
        #batch: vector or factor with batch label of each cell
        pca.data <- prcomp(xx[,-c(1:4)], center=FALSE) #compute PCA representation of the data
        yy <- pcRegression(pca.data,batch = xx$batch,n_top = top_pc)
        save(pca.data,yy,xx,top_pc,file = "pca_batch_reg.rda")
        return(list(pcr=yy,name=y@ID))
        return(yy)
    })

    save(res,file="res.rda")

    fres <- lapply(res, function(a){
        r2 <- as.data.frame(a$pcr$r2)
        r2 <- r2[1:min(top_pc,nrow(r2)),]
        #r2$R.squared <- cell_spec(r2$R.squared, bold = T,
        #                          color = ifelse(r2$p.value.lm <= 0.05, "red", "black"),
        #                          font_size = spec_font_size(r2$R.squared))
        #r2$R.squared <- cell_spec(r2$R.squared, )
        r2$PC <- row.names(r2)
        r2$PC <- as.integer(str_replace_all(r2$PC,pattern = "PC",replacement = ""))
        r2$dataSet <- a$name
        r2$ExplainedVar <- a$pcr$ExplainedVar[1:nrow(r2)]
        return(r2)
    })

    fres <- bind_rows(fres)
    fres$R.squared <- sprintf("%.3f",fres$R.square) %>% as.numeric()
    fres$R.squared <- cell_spec(fres$R.squared, bold = T,
                              color = ifelse(fres$p.value.lm <= 0.05, "red", "black"),
                              font_size = spec_font_size(fres$R.squared))
    fres2 <- fres %>% select(PC,R.squared,dataSet) %>% spread(key=dataSet,value=R.squared)

    ## explained_var
    explained_var <- fres %>% select(PC,ExplainedVar,dataSet) %>% spread(key=dataSet,value=ExplainedVar)

    ## pcRegscale
    pcRegscale <- sapply(res, function(a){
        return(a$pcr$pcRegscale)
    })
    pcRegscale <- data.frame(dataSet=names(pcRegscale),
                             pcRegscale=pcRegscale,
                             stringsAsFactors = FALSE)
    row.names(pcRegscale) <- NULL

    return(list(pcr=res,
                table=fres2,
                explained_var=explained_var,
                pcRegscale=pcRegscale,
                top_pc=top_pc))
}

get_pcr_table=function(x,top_pc=10){
    fres <- lapply(x, function(a){
        r2 <- as.data.frame(a$pcr$r2)
        r2 <- r2[1:min(top_pc,nrow(r2)),]
        #r2$R.squared <- cell_spec(r2$R.squared, bold = T,
        #                          color = ifelse(r2$p.value.lm <= 0.05, "red", "black"),
        #                          font_size = spec_font_size(r2$R.squared))
        #r2$R.squared <- cell_spec(r2$R.squared, )
        r2$PC <- row.names(r2)
        r2$dataSet <- a$name
        return(r2)
    })

    fres <- bind_rows(fres)
    return(fres)
}


calc_batch_effect_metrics=function(x,out_dir="./",prefix="test",missing_value_cutoff=0.5){

    cat("Missing value imputation ...\n")
    x <- lapply(x, function(y){
        y <- metaX::filterPeaks(y,ratio = missing_value_cutoff)
        y <- missingValueImpute(y)
        return(y)
    })
    #kbet_res <- run_kbet(x)
    sil_res <- calc_silhouette_width(x)
    pcr_res <- calc_pca_batch_regression(x)
    #res <- list(kbet=kbet_res,sil=sil_res,pcr=pcr_res)
    res <- list(kbet=NULL,sil=sil_res,pcr=pcr_res)
    return(res)
}


plot_pca=function(x,out_dir="./",prefix="test",pointSize=0.8,labelSize=4,
                  legendRowBatch = 10,
                  legendRowClass = NULL,show_class=FALSE,
                  pc="12",ncol_sub_fig=3,fig_res=120, pngWidth=6.5, legend_pos=NULL){
    # pc = "12" or "13"

    plotData <- lapply(x,function(y){y$pca_plot_data}) %>% bind_rows()

    if(pc == "12"){
        ggobj <- ggplot(data = plotData,aes(x=x,y=y,color=batch))
        ggobj <- ggobj + geom_hline(yintercept=0,colour="gray")+
            geom_vline(xintercept=0,colour="gray")+
            #geom_point()+
            #xlab(paste("PC1"," (",sprintf("%.2f%%",100*pca.res@R2[1]),") ",sep=""))+
            #ylab(paste("PC2"," (",sprintf("%.2f%%",100*pca.res@R2[2]),") ",sep=""))+
            xlab("PC1")+
            ylab("PC2")+
            theme_bw()+
            theme(#legend.position = legendPosition,
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
    }else{
        ggobj <- ggplot(data = plotData,aes(x=x,y=z,color=batch))
        ggobj <- ggobj + geom_hline(yintercept=0,colour="gray")+
            geom_vline(xintercept=0,colour="gray")+
            #geom_point()+
            #xlab(paste("PC1"," (",sprintf("%.2f%%",100*pca.res@R2[1]),") ",sep=""))+
            #ylab(paste("PC2"," (",sprintf("%.2f%%",100*pca.res@R2[2]),") ",sep=""))+
            xlab("PC1")+
            ylab("PC3")+
            theme_bw()+
            theme(#legend.position = legendPosition,
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
    }
    ggobj <- ggobj + stat_ellipse(aes(),geom = "path")
    if(show_class==FALSE){
        ggobj <- ggobj + geom_point(size=pointSize)
    }else{
        ggobj <- ggobj + geom_point(aes(shape=class),size=pointSize)+
            scale_shape_manual(values=1:n_distinct(plotData$class))
    }

    ggobj <- ggobj + facet_wrap(.~dataSet,ncol = ncol_sub_fig,scales="free")

    if(!is.null(legendRowBatch)){
        ggobj <- ggobj + guides(color = guide_legend(nrow = legendRowBatch))

    }

    if(!is.null(legend_pos)){
        ggobj <- ggobj + theme(legend.position=legend_pos)
    }

    #pdf(file = highfig,width = pdfWidth,height = pdfHeight)
    #pngWidth <- 6.5
    ht <- ceiling(length(unique(plotData$dataSet))/ncol_sub_fig) * pngWidth/ncol_sub_fig
    fig <- paste(out_dir,"/",prefix,"-pca_batch.png",sep="")
    png(fig,width = pngWidth,height = ht,res = fig_res,units = "in")

    print(ggobj)
    dev.off()

    return(list(fig=fig,data=plotData))

}


plot_density=function(x,out_dir="./",prefix="test",filter_by_quantile=0){
    a <- lapply(x, function(y){
        dat <- y@peaksData %>% mutate(dataSet=y@ID,sample=as.character(sample)) %>%
            as_tibble()
        dat$value[dat$value<=0] <- NA
        return(dat)
    })

    if(filter_by_quantile>0){
        a <- lapply(a, function(y){
            b <- quantile(y$value,probs = c(filter_by_quantile,1-filter_by_quantile),na.rm = TRUE)
            y <- filter(y,value >= b[1], value <= b[2])
            return(y)
        })
    }

    dat <- bind_rows(a)
    #dat$value[dat$value<=0] <- NA

    gg <- dat %>% ggplot(aes(x=log2(value),color=sample)) +
        #stat_density(geom="line")+
        geom_line(stat="density",alpha=0.4)+
        #geom_density()+
        facet_wrap(.~dataSet,ncol = 3,scales="free")+
        guides(color=FALSE)+
        theme_bw()+
        theme(#legend.position = legendPosition,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

    pngWidth <- 6.5
    ht <- ceiling(length(unique(dat$dataSet))/3) * pngWidth/3
    fig <- paste(out_dir,"/",prefix,"-density.png",sep="")
    png(fig,width = pngWidth,height = ht,res = 120,units = "in")
    print(gg)
    dev.off()
    return(fig)
}


calc_ml_metrics=function(x,sample_list=NULL,sample_class=NULL,use_all=TRUE,
                         missing_value_cutoff=0.5, valueID="value",cpu=0,
                         n_repeats=20){

    if(missing_value_cutoff > 0){
        x <- lapply(x, filterPeaks, ratio=missing_value_cutoff)
    }

    if(use_all==FALSE){
        ## use common genes/proteins
        for(i in 1:length(x)){
            if(i == 1){
                common_ids <- x[[1]]@peaksData$ID
            }else{
                common_ids <- intersect(common_ids,x[[i]]@peaksData$ID)
            }
        }
        common_ids <- unique(common_ids)
        cat("Common IDs:", length(common_ids),"\n")

        x <- lapply(x, function(y){
            y@peaksData <- y@peaksData %>% filter(ID %in% common_ids)
            return(y)
        })

    }

    if(!is.null(sample_list)){
        sam <- read.delim(sample_list,stringsAsFactors = FALSE) %>%
            dplyr::select(sample,class)
        x <- lapply(x, function(y){
            y@sampleList <- read.delim(sample_list,stringsAsFactors = FALSE)
            y@peaksData$class<- NULL
            y@peaksData <- merge(y@peaksData,sam,by="sample")
            return(y)
        })

        class_group <- sam$class %>% unique()
    }else{
        if(!is.null(sample_class)){
            class_group <- sample_class
        }else{
            class_group <- x[[1]]@peaksData$class %>% unique()
        }
    }

    cat("Use class:",paste(class_group,collapse = ","),"\n")

    res <- lapply(x, function(y){

        if(requireNamespace("doMC",quietly = TRUE)){
            if(cpu==0){
                library(doMC)
                cpu = detectCores()
                registerDoMC(cores = cpu)
            }else if(cpu>=2){
                library(doMC)
                registerDoMC(cores = cpu)
            }
        }
        if(any(is.na(y@peaksData$value))){
            y <- metaX::missingValueImpute(y)
        }else{
            cat("No missing value!\n")
        }


        dd <- lapply(1:n_repeats,function(i){
            cat("Repeat:", i,"\n")
            dat <- featureSelection(y,group=class_group,method = "rf",
                                    valueID = valueID, fold = 5,
                                    resampling_method = "LOOCV",
                                    repeats = 10,
                                    plot_roc = FALSE,
                                    ratio = 2/3, k = 100,
                                    metric = "ROC",
                                    sizes = length(unique(y@peaksData$ID)),
                                    plotCICurve = FALSE,verbose=TRUE)
            dat$results$Repeat_ID = i
            return(dat)

        })

        eval_res <- lapply(dd,function(x){ x$results %>% select(everything())}) %>%
            bind_rows()

        return(eval_res)
        #return(dat)
    })

    name_datasets <- names(res)

    for(i in 1:length(res)){
        res[[i]]$dataSet <- name_datasets[i]
    }

    res_table <- bind_rows(res)

    ftable <- res_table %>% group_by(dataSet) %>% dplyr::summarise(mean_ROC=mean(ROC),median_ROC=median(ROC),sd_ROC=sd(ROC))

    fres <- list(data=res_table,table=ftable,class_group=class_group)

    return(fres)
}

plot_ml_boxplot=function(x, out_dir="./",prefix="test"){

    mean_roc <- x %>% group_by(dataSet) %>% summarise(mean_ROC=mean(ROC)) %>% arrange(desc(mean_ROC))
    x$dataSet <- factor(x$dataSet,levels = mean_roc$dataSet)



    fig <- paste(out_dir,"/",prefix,"-ml_boxplot.png",sep="")
    png(fig,width = 800,height = 400,res=120)
    gg <- ggplot(x,aes(x=dataSet,y=ROC))+
        geom_boxplot(width=0.5,outlier.size=0.2)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        stat_summary(fun.y=mean, colour="darkred", geom="point",
                     shape=18, size=3, show.legend = FALSE) +
        geom_text(data = mean_roc %>% mutate(mean_roc_label=sprintf("%.2f",mean_ROC)),
                  aes(label = mean_roc_label, y = mean_ROC + 0.1*(max(x$ROC)-min(x$ROC))),colour="darkred")
    print(gg)
    dev.off()

    ## save plot data
    fig_data_file <- paste(out_dir,"/",prefix,"-ml_boxplot.rds",sep="")
    fig_data <- list()
    fig_data$plot_object <- gg
    fig_data$plot_data <- x
    saveRDS(fig_data,file = fig_data_file)

    return(fig)

}

# generate overview table and overview radar figure
generate_overview_table=function(x,highlight_top_n=3,min_auc=0.8){
    ## generate an overview table which contains different metrics
    ## input: x =>

    dat <- get_identification_summary_table(x,format = FALSE)

    if(!is.null(x$batch_effect_metrics)){

        # kBET
        if(!is.null(x$batch_effect_metrics$kbet)){
            dat <- merge(dat,x$batch_effect_metrics$kbet$table %>%
                             select(dataSet,kBET.observed))

            dat$kBET <- 1 - dat$kBET.observed
            dat$kBET.observed <- NULL
        }
        #dat$kBET.observed <- cell_spec(dat$kBET.observed,
        #                               color = ifelse(y >= max(y), "red", "black"))


        # Silhouette width
        sil <- data.frame(dataSet=names(x$batch_effect_metrics$sil),
                          silhouette_width=x$batch_effect_metrics$sil)
        dat <- merge(dat,sil)
        dat$silhouette_width <- abs(dat$silhouette_width)

        # PCR
        # pcr <- get_pcr_table(x$batch_effect_metrics$pcr$pcr)
        # sig_pc <- pcr %>% dplyr::filter(p.value.lm<=0.05) %>%
        #     dplyr::select(PC) %>%
        #     dplyr::distinct()
        # pcr <- pcr %>% dplyr::filter(PC %in% sig_pc$PC) %>%
        #     dplyr::select(dataSet,PC,R.squared) %>%
        #     dplyr::mutate(PC=paste("batch_effect_",PC,sep="")) %>%
        #     tidyr::spread(key = PC,value = R.squared)
        pcRegscale <- x$batch_effect_metrics$pcr$pcRegscale
        pcRegscale$pcRegscale <- 1 - pcRegscale$pcRegscale

        dat <- merge(dat,pcRegscale)
    }

    if(!is.null(x$network_table)){
        dat <- merge(dat,x$network_table$cor %>% dplyr::select(dataSet,ks) %>%
                         dplyr::rename(complex_ks=ks))
    }




    if(!is.null(x$input_parameters$x2)){
        dat <- merge(dat,x$protein_rna$feature_wise_cor_table %>%
                         dplyr::select(dataSet,median_cor) %>%
                         dplyr::rename(gene_wise_cor=median_cor))

        dat <- merge(dat,x$protein_rna$sample_wise_cor_table %>%
                         dplyr::select(dataSet,median_cor) %>%
                         dplyr::rename(sample_wise_cor=median_cor))
    }

    if(!is.null(x$ml)){
        dat <- merge(dat,x$ml$table %>%
                         dplyr::select(dataSet,mean_ROC) %>%
                         dplyr::rename(AUROC=mean_ROC))
    }

    if(!is.null(x$fun_pred)){
        f_res <- get_func_pred_meanAUC(x$fun_pred$data,min_auc = min_auc)
        dat <- merge(dat,f_res)
    }

    return(dat)

}



plot_radar=function(x){

    radar_dat <- x

    p <- plot_ly(
        type = 'scatterpolargl',
        #fill = 'toself',
        #alpha = 0.1,
        #line = list(smoothing = 1, shape = "spline"),
        mode = 'lines+markers'#,
        #marker = list(color = "white",
        #              size = 4,
        #              line = list(width = 2)),

    )
    #
    for(i in 1:nrow(radar_dat)){
        p <- p %>% add_trace(
            r = radar_dat[i,-1] %>% as.numeric(),
            theta = names(radar_dat)[-1],
            name = x$dataSet[i]
            #line = list(color = "#d3d3d3", dash = "3px")
            #mode = "lines"


        )
    }
    p <- p %>%
        layout(
            polar = list(
                radialaxis = list(
                    visible = T,
                    range = c(0,1)
                )
            )
        )

    return(p)
}




