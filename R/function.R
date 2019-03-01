
calc_function_prediction_metrics=function(x,missing_value_cutoff=0.5,
                                          use_all=TRUE,
                                          sample_list=NULL,
                                          sample_class=NULL,
                                          cpu=0,
                                          out_dir="./",
                                          prefix="test",
                                          min_auc=0.8){
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



    if(cpu==0){
        cpu <- detectCores()
    }
    input_cpu <- cpu
    if(cpu > length(x)){
        cpu <- length(x)
    }
    cat("Used cpus:",cpu,"\n")
    cl <- makeCluster(getOption("cl.cores", cpu))
    clusterExport(cl, c("MatNet"),envir=environment())
    clusterExport(cl, c("getTable"),envir=environment())
    #save(x,file="x.rda")

    net_res <- parLapply(cl,x,function(y){
        y@peaksData <- y@peaksData %>% dplyr::filter(class %in% class_group)
        dat <- getTable(y,valueID="value")
        row.names(dat) <- dat$ID
        dat$ID <- NULL
        dat <- as.matrix(dat)
        net_data <- MatNet(dat)
        return(net_data)
    })

    stopCluster(cl)
    #save(net_res,file="net_res.rda")


    fun_res <- lapply(net_res, function(y){
        pres <- function_predict(y,cpu=input_cpu)
    })

    dat_name <- names(fun_res)
    for(i in 1:length(dat_name)){
        fun_res[[i]]$dataSet <- dat_name[i]

    }

    fres <- bind_rows(fun_res)

    f_table <- fres
    f_table$AUC[is.na(f_table$AUC)] <- 0

    f_table$AUC <- sprintf("%.3f",f_table$AUC) %>% as.numeric()
    # max4term <- f_table[,-c(1,2)] %>% apply(1, max,na.rm = TRUE)
    #f_table <- f_table[max4term>=min_auc,]

    f_table_format <- f_table %>% dplyr::select(dataSet,term,AUC,db_num) %>%
        spread(key=dataSet,value=AUC)

    max4term <- f_table_format[,-c(1,2)] %>% apply(1,max)
    f_table_format <- f_table_format[max4term>=min_auc,]

    ## generate boxplot
    fig <- plot_fun_boxplot(f_table_format[,-c(1,2)],out_dir = out_dir,
                            prefix = prefix)

    f_table_format2 <- f_table_format %>% mutate(db_num=NULL)
    row.names(f_table_format2) <- f_table_format2$term
    f_table_format2$term <- NULL
    f_table_format2 <- t(f_table_format2)

    final_table <- f_table_format2

    for(i in 1:ncol(f_table_format2)){
        y <- f_table_format2[,i]
        y <- cell_spec(y, bold = T,
                       color = ifelse(y >= max(y), "red", "black"),
                       font_size = spec_font_size(y,scale_from=range(f_table$AUC)))
        final_table[,i] <- y

    }

    final_table <- t(final_table)

    return(list(data=fres,table=final_table,fig=fig))

}

function_predict=function(net_data,min_n=10,kfold=5,ranNum=1000,r=0.5,cpu=0){

    require(igraph)
    require(pROC)
    net_data <- as.matrix(net_data)
    net_data <- graph.edgelist(net_data,directed=F)

    #goauc <- data.frame(go=flist,gonum=0,rnanum=0,rnaauc=0,stringsAsFactors=F)

    #goaocList <- list()

    fundb <- readRDS(system.file("extdata/kegg.rds",package = "OmicsEV"))
    #fundb <- fundb %>% dplyr::filter(n>=min_n)

    fundb_list <- list()
    for(i in 1:nrow(fundb)){
        fundb_list[[i]] <- fundb$ID[i]
    }
    names(fundb_list) <- fundb$fun

    if(cpu==0){
        cpu <- detectCores()
    }
    if(cpu > length(fundb_list)){
        cpu <- length(fundb_list)
    }
    cat("Used cpus:",cpu,"\n")
    cl <- makeCluster(getOption("cl.cores", cpu))
    clusterExport(cl, c(".kfoldCrossValidation"),envir=environment())
    clusterExport(cl, c("V"),envir=environment())
    clusterExport(cl, c("get.adjacency"),envir=environment())
    clusterExport(cl, c("roc"),envir=environment())
    clusterExport(cl, c("min_n"),envir=environment())
    #clusterExport(cl, c("net_data"),envir=environment())

    #save(fundb_list,net_data,file="test.rda")
    pres <- parLapply(cl,fundb_list,function(y){
        ann <- strsplit(x = y,split = ";") %>% unlist()

        n_ID1 <- length(ann)
        ann_use <- intersect(ann,V(net_data)$name)
        n_ID2 <- length(ann_use)

        if(n_ID2 >= min_n){

            res_auc <- tryCatch({
                res_auc <- .kfoldCrossValidation(geneset=ann_use,network=net_data,r=r,
                                                 kfold=kfold,ranNum=ranNum)
                res_auc
            },error=function(e){
                message("Please see the file: kfoldCrossValidation.rda for related data!")
                save(e,ann_use,net_data,r,kfold,ranNum,file="kfoldCrossValidation.rda")
                stop("error in .kfoldCrossValidation!")
                return(NULL)
            },
            warning=function(cond){
                message("Please see the file: kfoldCrossValidation.rda for related data!")
                save(e,ann_use,net_data,r,kfold,ranNum,file="kfoldCrossValidation.rda")
                return(NULL)
            })

            aucR <- c()

            for(j in c(1:kfold)){
                aucR <- c(aucR,res_auc[[j]]$auc)
            }

            fauc <- mean(aucR)

            dd <- data.frame(db_num=n_ID1,ID_num=n_ID2,AUC=fauc,ID=y,
                             stringsAsFactors = FALSE)
        }else{
            dd <- data.frame(db_num=n_ID1,ID_num=n_ID2,AUC=NA,ID=y,
                             stringsAsFactors = FALSE)
        }
        return(dd)

    })
    stopCluster(cl)

    dat_name <- names(pres)
    for(i in 1:length(dat_name)){
        pres[[i]]$term <- dat_name[i]

    }

    fres <- bind_rows(pres)

    return(fres)

}

.kfoldCrossValidation <- function(geneset,network,r=0.5,kfold=5,ranNum=1000){
    #save(geneset,network,r,kfold,ranNum,file = "kfold.rda")
    geneset <- sample(geneset,length(geneset))   ##random gene set

    ##divide geneset into kfold sets####
    m <- floor(length(geneset)/kfold)
    n <- m*kfold
    d <- length(geneset)-n
    label <- c()
    if(d!=0){
        for(i in c(1:kfold)){
            if(i<=d){
                label <- c(label,rep(i,m+1))
            }else{
                label <- c(label,rep(i,m))
            }
        }
    }else{
        for(i in c(1:kfold)){
            label <- c(label,rep(i,m))
        }
    }
    aur <- c()

    negativeSet <- setdiff(V(network)$name,geneset)

    foldList <- list()

    for(i in c(1:kfold)){

        validateSet <- geneset[which(label==i)]
        trainSet <- geneset[which(label %in% setdiff(c(1:kfold),i))]



        pt1 <- .netwalker(trainSet,network,r)
        gS <- data.frame(name=V(network)$name,label=0,score=pt1,globalP=(rank(-pt1,ties.method="min")-1)/length(pt1),stringsAsFactors=F)
        #gS <- data.frame(name=V(network)$name,label=0,score=pt1,globalP=(rank(-pt1,ties.method="min")-1)/length(pt1),localP=1,edgeP=1,maxP=1,stringsAsFactors=F)

        #calculate local P and edge P. For the local P, we first randomly select nodes with the same number of seeds and used these random nodes as seed to calcualte  the p values of all nodes.
        #For edge P, we first genenerate the random network with the same number of nodes and the same degree distribution and then use the same seeds to calculate the p values of all nodes based on the random network

        gS[gS[,1] %in% validateSet,2] <- 1
        gS[gS[,1] %in% negativeSet,2] <- 0

        #re <- .netwalker_random(trainSet,network,r,ranNum)

        #ranPt <- re$randomSeed
        #ranEg <- re$randomEdge

        #x <- cbind(pt1,ranPt)

        #localP <- apply(x,1,.calculateLocalP)
        #gS[,5] <- localP

        #x <- cbind(pt1,ranEg)
        #edgeP <- apply(x,1,.calculateLocalP)
        #gS[,6] <- edgeP
        #gS[,7] <- apply(gS[,c(4:6)],1,max)

        gS <- gS[gS[,1] %in% setdiff(gS[,1],trainSet),]

        roc1 <- roc(gS[,2],gS[,4])

        sensitivity <- roc1$sensitivities
        minusSpecificity <- 1-roc1$specificities

        p5 <- wilcox.test(sensitivity,minusSpecificity)
        p5 <- p5$p.value
        auc <- as.numeric(roc1$auc)

        re <- list(gS=gS,roc=roc1,auc=auc,p5=p5)

        foldList[[i]] <- re

        cat("Fold:",i,"\n",sep="")
    }

    return(foldList)
}

.calculateLocalP <- function(v){
    rP <- v[1]
    localP <- v[2:length(v)]
    p <- length(which(rP<localP))/length(localP)
    return(p)
}

.netwalker=function(seed, network, r=0.5){
    adjMatrix <- get.adjacency(network)
    adjMatrix <- as.matrix(adjMatrix)
    de <- apply(adjMatrix,2,sum)
    W <- t(t(adjMatrix)/de)

    p0 <- array(0,dim=c(nrow(W),1))
    rownames(p0) <- rownames(W)
    p0[seed,1] <- 1/length(seed)

    pt <- p0
    pt1 <- (1-r)*(W%*%pt)+r*p0

    while(sum(abs(pt1-pt))>1e-6){
        pt <- pt1
        pt1 <- (1-r)*(W%*%pt)+r*p0
    }

    return(pt1)
}

.netwalker_random <- function(seed,network,r=0.5,ranNum=1000){

    ##local P:  the random walker start from the same number of randomly selected transcription factors to calculate random scores for each node

    node <- V(network)$name
    degree <- igraph::degree(network)
    randomRe <- array(0,dim=c(length(node),ranNum))
    randomEdge <- array(0,dim=c(length(node),ranNum))

    for(i in c(1:ranNum)){
        if(i %% 100 == 0){
            cat("Random:",i,"\n",sep="")
        }

        ###random seed

        ranSeed <- sample(node,length(seed))
        ranScore <- .netwalker(ranSeed,network,r=r)
        randomRe[,i] <- ranScore

        ####random edge

        ranNetwork <- degree.sequence.game(degree,method="vl")
        V(ranNetwork)$name <- node
        ranScore <- .netwalker(seed,ranNetwork,r=r)
        randomEdge[,i] <- ranScore

    }
    re <- list(randomSeed=randomRe,randomEdge=randomEdge)
    return(re)
}


plot_fun_boxplot=function(x, out_dir="./",prefix="test"){

    rank_x <- apply(x, 1, function(y){rank(-y,ties.method = "min")}) %>% t
    ## median or mean
    sort_name <- names(sort(apply(rank_x, 2, mean)))
    x <- x[,sort_name]
    rank_x <- rank_x[,sort_name]

    box_data <- gather(as.data.frame(rank_x),"Method","Rank")
    box_data$Method <- factor(box_data$Method,levels = unique(box_data$Method))

    fig <- paste(out_dir,"/",prefix,"-fun_boxplot.png",sep="")
    png(fig,width = 800,height = 400,res=120)
    gg <- ggplot(box_data,aes(x=Method,y=Rank))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        stat_summary(fun.y=mean, colour="darkred", geom="point",
                     shape=18, size=3, show.legend = FALSE) +
        geom_text(data = box_data %>% group_by(Method) %>%
                      dplyr::summarise(mean_rank=mean(Rank)) %>%
                      mutate(mean_rank_label=sprintf("%.3f",mean_rank) ),
                  aes(label = mean_rank_label, y = mean_rank + 0.5),colour="darkred")
    print(gg)
    dev.off()
    return(fig)

}

get_func_pred_table=function(x, min_auc=0.8){
    f_table <- x
    f_table$AUC[is.na(f_table$AUC)] <- 0

    f_table$AUC <- sprintf("%.3f",f_table$AUC) %>% as.numeric()
    # max4term <- f_table[,-c(1,2)] %>% apply(1, max,na.rm = TRUE)
    #f_table <- f_table[max4term>=min_auc,]

    f_table_format <- f_table %>% dplyr::select(dataSet,term,AUC,db_num) %>%
        tidyr::spread(key=dataSet,value=AUC)

    max4term <- f_table_format[,-c(1,2)] %>% apply(1,max)
    f_table_format <- f_table_format[max4term>=min_auc,]

    auc_data <- f_table_format[,-c(1,2)]
    rank_x <- apply(auc_data, 1, function(y){rank(-y,ties.method = "min")}) %>% t
    ## median or mean
    sort_name <- names(sort(apply(rank_x, 2, mean)))
    auc_data <- auc_data[,sort_name]
    rank_x <- rank_x[,sort_name]

    box_data <- tidyr::gather(as.data.frame(rank_x),"Method","Rank")
    box_data$Method <- factor(box_data$Method,levels = unique(box_data$Method))

    res <- box_data %>% group_by(Method) %>%
        summarise(fun_rank=mean(Rank)) %>%
        rename(dataSet=Method) %>% ungroup()
    return(res)
}


do_function_pred=function(data_dir,sample_list=NULL,
                          missing_value_cutoff=0.5,
                          use_common_features_for_func_pred=FALSE,
                          cpu=0,
                          out_dir="./",
                          prefix="omicsev"){
    ## import data
    input_data_files <- list.files(path = data_dir,pattern = ".tsv",
                                   full.names = TRUE,include.dirs = TRUE)

    ## the number of datasets
    dataset_names <- basename(input_data_files) %>%
        str_replace_all(pattern = ".tsv",replacement = "")

    x1 <- list()
    for(i in 1:length(input_data_files)){
        x1[[i]] <- import_data(input_data_files[i],sample_list = sample_list,
                               data_type="gene",
                               missing_value_cutoff = missing_value_cutoff)
        x1[[i]]@ID <- dataset_names[i]
    }
    names(x1) <- dataset_names
    fp_res <- calc_function_prediction_metrics(x1,
                                               missing_value_cutoff=missing_value_cutoff,
                                               use_all=!use_common_features_for_func_pred,
                                               cpu=cpu,
                                               out_dir=out_dir,
                                               prefix=prefix)
    saveRDS(fp_res,file = paste(out_dir,"/res.rds",sep=""))
    return(fp_res)
}




