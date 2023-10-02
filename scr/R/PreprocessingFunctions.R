library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(paletteer)

# ========== define function ==========
## === DEGs ===
obtain_DEGs <- function(rawdata_path, group_path, 
                        rpkm_path, thre = 1,
                        output_path.rawDEGs, DEGs.filtered.path){
  "
  :para rawdata_path: pathway to rawdata, it could be OTU tables, reads number of metagenomics and metatranscriptomics, etc.
  :para group_path: pathway to group information, need be create manually by user
      NOTE! When making group tables, please make sure the control is placed on the first place. 
  :para output_path: pathway to save all DESeq results
  "
  # load rawdata
  raw.counts <- read.table(rawdata_path, sep='\t', header = T,row.names = 1)
  
  # load group information
  group <- read.csv(group_path)
  
  # --- read file directly | DESeq analysis ---
  file_path <- output_path.rawDEGs
  if (file.exists(file_path)) {
    print(paste("File", file_path, "exists."))
    DEGs_raw <- read.csv(file_path)
  } else {
    # Differential analysis
    set.seed(0)
    suppressPackageStartupMessages(library(DESeq2))
    dds_raw <- DESeqDataSetFromMatrix(raw.counts, group, design = ~Group)
    dds <- DESeq(dds_raw)
    
    # set control and obtain info. of DEGs
    conditions <- unique(group$Group)
    control <- conditions[1] 
    loc_control <- match(control, conditions)
    conditions <- conditions[-loc_control]
    
    # merge DEGs reults of all groups
    for (ref in conditions){
      DESeqRes <- results(dds,contrast = c("Group", ref, control)) 
      DEGs <- as.data.frame(DESeqRes)
      DEGs$degs <- paste(ref, 'vs', control, sep = '_')
      DEGs$id <- rownames(DEGs)
      # remove row with NA
      DEGs <- na.omit(DEGs)
      if(ref == conditions[1]){
        DEGs_raw <- DEGs
      }else{
        DEGs_raw <- rbind(DEGs_raw, DEGs)
      }
    }
    # --- save data: raw DEGs ---
    write.csv(DEGs_raw, output_path.rawDEGs, row.names = FALSE)
  }
  
  # --- obtain filtered DEGs ---
  file_path <- DEGs.filtered.path
  if (file.exists(file_path)) {
    print(paste("File", file_path, "exists."))
    DEGs.filtered.merged <- read.csv(file_path)
  } else {
    # === Filter Differential expressed genes (DEGs) ===
    # --- threshold 1: padj < 0.05, |log2 fold change| > 1 --- 
    print(paste('The total number of DEGs (not unique):',length(DEGs_raw$id)))
    print(paste('The total number of unique DEGs (not unique):', length(unique(DEGs_raw$id))))
    # DEGs.filtered <- subset(DEGs_raw,  padj < 0.05 & abs(log2FoldChange) > 1) ## padj
    DEGs.filtered <- subset(DEGs_raw,  pvalue < 0.05 & abs(log2FoldChange) > 1) ## padj
    print(paste('The total number of deseq filtered DEGs:',length(DEGs.filtered$id))) 
    print(paste('The total number of unique deseq filtered DEGs:', length(unique(DEGs.filtered$id)))) 
    
    # --- threshold 2: mean > 1/length ---
    # obtian expression levels (fpkm)
    rpkm <- read.table(rpkm_path, sep='\t', header = T, row.names = 1)
    rpkm$mean <- apply(rpkm, 1, mean)
    rpkm$id <- rownames(rpkm)
    print(paste('The number of total genes (unique):',length(unique(rpkm$id))))
    loc <- rpkm$mean >= thre # 1/(length(rpkm)-2)
    rpkm.filtered <- rpkm[loc,]
    
    # --- merge filtering conditions ---
    DEGs.filtered.merged <- merge(DEGs.filtered, rpkm.filtered, by = 'id', all = F)
    print(paste('The number of total DEGs after filtering :',length(DEGs.filtered.merged$id)))
    print(table(DEGs.filtered.merged$degs))
    
    # # merge filtering conditions
    # df.filtered <- merge(DEGs.filtered, rpkm.filtered, by = 'id', all = F)
    # DEGs_BvsD <- subset(df.filtered, df.filtered$degs == 'Blue_vs_Dark')
    # DEGs_BvsD$B_DEGs_label <- 1
    # DEGs_BvsD$Y_DEGs_label <- 0
    # 
    # DEGs_YvsD <- subset(df.filtered, df.filtered$degs == 'Yellow_vs_Dark')
    # DEGs_YvsD$B_DEGs_label <- 0
    # DEGs_YvsD$Y_DEGs_label <- 1
    # 
    # DEGs.filtered.merged <- merge(DEGs_BvsD, DEGs_YvsD, by = 'id', all = T)[,c(1:(8+length(rpkm)+1))]
    # colnames(DEGs.filtered.merged) <- colnames(DEGs_BvsD)
    # print(paste('The number of total DEGs after filtering :',length(DEGs.filtered.merged$id)))
    # --- save data: filtered DEGs ---
    write.csv(DEGs.filtered.merged, DEGs.filtered.path, row.names = FALSE)
  }
  return(DEGs.filtered.merged)
}


get_tsne_df <- function(DEGs.filtered, dataset,
                        group_path, add_info=FALSE,
                        anno.info.1=FALSE, anno.info.2=FALSE,
                        perplex = 20){ # anno.info.1=anno.signal # anno.info.2=anno.tmhmm
  # === data for tsne ===
  # head(DEGs.filtered) # checking code
  df <- unique(DEGs.filtered[,-c(2:7, length(DEGs.filtered))])
  # head(df) # checking code
  
  # normalize the expression levels 
  df.norm <- log(df[,-(1:2)] + 1, 10)
  df.norm$id <- df$id
  df.norm$degs <- df$degs
  
  # --- grouping to different datasets ---
  DEG_group <- as.vector(as.data.frame(table(DEGs.filtered$degs))$Var1)
  
  # load conditions information
  # set control and obtain info. of DEGs
  conditions <- read.csv(group_path)
  # conditions <- unique(sample_group$Group)
  # control <- conditions[1] 
  # loc_control <- match(control, conditions)
  # conditions <- conditions[-loc_control]
  
  group <- DEG_group[1]
  mylist <- list()
  file_path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
  if (file.exists(file_path)) {
    for (i in 1:length(DEG_group)){
      group <- DEG_group[i]
      file_path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
      graph_data <- read.csv(file_path)
      mylist[group] <- list(graph_data)
      print(paste("File", file_path, "exists and successfully load."))
    }
  } else {
    for (i in 1:length(DEG_group)){
      group <- DEG_group[i]
      tsne.temp <- subset(df.norm, df.norm$degs == group) # unique genes, but some might share the same expresion pattern, such as 0*n. tsne.temp: expression, degs label, id
      
      ### ==TODO: 提前unique到targeted samples== ###
      tsne.raw <- tsne.temp[,1:(length(tsne.temp)-2)]
      print(paste('The total number of unique DEGs for tsne analysis of', group, 'dataset is:', length(tsne.temp$id)))
      
      condition <- as.vector(
        as.data.frame(table(conditions$Group[-(1:3)]))$Var1[i]
      )
      samples <- subset(conditions, Group == condition | Group == conditions$Group[1])$Samples
      df_tsne <- tsne.raw[,samples]
      # obtain unique expression pattern
      df_tsne.unique <- unique(df_tsne)
      # obtain none-zero genes
      tsne.conditions <- tsne.temp[c(samples, 'id', 'degs')]
      tsne.No0 <- tsne.conditions[as.logical(rowSums(tsne.conditions[,1:(length(tsne.conditions)-1)] != 0)), ]
      
      # --- tsne analysis ---
      library(Rtsne)
      set.seed(0)
      tsne_fpkm = Rtsne(
        tsne.No0[,1:6],
        dims = 2,
        pca = T,
        max_iter = 1000,
        theta = 0.4,
        perplexity = perplex,
        verbose = F
      )
      # merge tsne results
      tsne_result = as.data.frame(tsne_fpkm$Y)
      colnames(tsne_result) = c("tSNE1","tSNE2")
      
      # # --- tsne check ---
      # ggplot(tsne_result,aes(tSNE1,tSNE2)) +  # ,colour = degs
      #   geom_point(alpha = 0.5, size = 1) + 
      #   # geom_point(colour="#85A4FE", size = 1) + 
      #   labs(title = '') + 
      #   theme_bw() + 
      #   scale_colour_manual(values = pale_25) +
      #   mytheme1
      
      # --- construct graph data ---
      # load additional information: Subcellular information including SignalP and tmhmm
      if (anno.info.1 != FALSE) {
        sub.signalp <- read.csv(anno.info.1, header = T)
        sub.tmhmm <- read.csv(anno.info.2, header = T)
        subloca <- merge(sub.signalp[, c('id', 'X.')], sub.tmhmm[, c('id', 'Topology')], by = 'id')
        ## signalp
        loc <- which(subloca$X. == 'N') 
        subloca$Signalp[loc] <- 0
        loc <- which(subloca$X. != 'N')
        subloca$Signalp[loc] <- 1
        table(subloca$Signalp)
        ## tmhmm
        loc <- which(subloca$Topology == 'Topology=o') 
        subloca$tmhmm[loc] <- 0
        loc <- which(subloca$Topology != 'Topology=o')
        subloca$tmhmm[loc] <- 1
        table(subloca$tmhmm)
      }
      
      # --- combine to obtain graph data (all info) & dgi data--- 
      # tsne_info <- cbind(tsne_result, df_tsne, tsne.temp[,c('degs', 'id')])
      tsne_info <- cbind(tsne_result, tsne.No0)
      
      if (anno.info.1 != FALSE) {
        graph_data <- merge(tsne_info, subloca, by.y = 'id')
        dgi_data <- graph_data[,c((2:(length(df_tsne)+2)), (length(graph_data)-1):length(graph_data))]
        rownames(dgi_data) <- graph_data$id
      } else {
        graph_data <- tsne_info
        dgi_data <- graph_data[,c((1:(length(df_tsne)+2)))]
        rownames(dgi_data) <- graph_data$id
      }
      
      # --- save data ---
      graphdata_path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
      dgidata_path <- paste("./data/", dataset, "/inputdata/dgidata_", group, '.csv', sep = '')
      write.csv(graph_data, graphdata_path, row.names = FALSE)
      write.csv(dgi_data, dgidata_path, row.names = FALSE)
      
      # temp variable
      mylist[group] <- list(graph_data)
    }
  }
  
  return(mylist)
}

# dimension reduction for more option
reduce_df <- function(DEGs.filtered, dataset,
                      group_path, add_info=FALSE,
                      anno.info.1=FALSE, anno.info.2=FALSE,
                      Reduction='tsne', # or 'umap'
                      normalization='log', # or 'stand'
                      perplex = 20){ # anno.info.1=anno.signal # anno.info.2=anno.tmhmm
  # --- extract data for dimension redcution ---
  df <- unique(DEGs.filtered[,-c(2:7, length(DEGs.filtered))])
  
  # --- normalization method ---
  if (Reduction == 'log'){
    # normalize the expression levels 
    df.norm <- log(df[,-(1:2)] + 1, 10)
    df.norm$id <- df$id
    df.norm$degs <- df$degs
  } else {
    # import data
    rawdata <- df[,3:length(df)]
    # Standardize the data
    df.norm <- data.frame(t(apply(rawdata, 1, function(v) {
      (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    })), stringsAsFactors = FALSE)
    # Replace NA values with 0
    df.norm[is.na(df.norm)] <- 0
    df.norm$id <- df$id
    df.norm$degs <- df$degs
  }
  
  # --- grouping to different datasets ---
  DEG_group <- as.vector(as.data.frame(table(DEGs.filtered$degs))$Var1)
  
  # --- load conditions information ---
  # set control and obtain info. of DEGs
  conditions <- read.csv(group_path)
  group <- DEG_group[1]
  
  # --- obtain graph_data list ---
  mylist <- list()
  file_path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
  if (file.exists(file_path)) {
    for (i in 1:length(DEG_group)){
      group <- DEG_group[i]
      file_path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
      graph_data <- read.csv(file_path)
      mylist[group] <- list(graph_data)
      print(paste("File", file_path, "exists and successfully load."))
    }
  } else {
    for (i in 1:length(DEG_group)){
      # --- obtain the first paried group ---
      group <- DEG_group[i]
      # obtain the samples and condition information of targeted groups
      condition <- as.vector(
        as.data.frame(table(conditions$Group[-(1:3)]))$Var1[i]
      )
      samples <- subset(conditions, Group == condition | Group ==conditions$Group[1])$Samples
      
      # --- obtain data of targeted group for dimension reduction ---
      df.temp <- subset(df.norm, df.norm$degs == group) # unique genes, but some might share the same expresion pattern, such as 0*n. df.temp: expression, degs label, id
      df.data <- df.temp[,1:(length(df.temp)-2)]
      print(paste('The total number of unique DEGs for tsne analysis of', group, 'dataset is:', length(df.temp$id)))
      df.data.group <- df.data[,samples]
      # obtain unique expression pattern
      df.data.group.unique <- unique(df.data.group)
      # obtain none-zero genes
      df.group.unique.info <- df.temp[c(samples, 'id', 'degs')]
      df.No0 <- df.group.unique.info[as.logical(rowSums(df.group.unique.info[,1:(length(df.group.unique.info)-2)] != 0)), ]
      
      # --- choose dimension reduction method ---
      if (Reduction == 'umap'){
        # UMAP analysis
        library(umap)
        set.seed(0)
        umap <- umap(df.No0[,1:6], method = 'naive', n_neighbors = 10)
        
        # Create a dataframe for plotting
        df.reduced <- data.frame(umap$layout)
        colnames(df.reduced) <- c('umap1', 'umap2')
        
      } else {
        # tsne analysis
        library(Rtsne)
        set.seed(0)
        tsne_fpkm = Rtsne(
          df.No0[,1:6],
          dims = 2,
          pca = T,
          max_iter = 1000,
          theta = 0.4,
          perplexity = perplex,
          verbose = F
        )
        
        # merge tsne results
        df.reduced <- as.data.frame(tsne_fpkm$Y)
        colnames(df.reduced) = c("tSNE1","tSNE2")
      }
      
      df.reduced.info <- cbind(df.reduced, df.No0)
      
      # --- subcellular information --
      # load additional information: Subcellular information including SignalP and tmhmm
      if (anno.info.1 != FALSE) {
        sub.signalp <- read.csv(anno.info.1, header = T)
        sub.tmhmm <- read.csv(anno.info.2, header = T)
        subloca <- merge(sub.signalp[, c('id', 'X.')], sub.tmhmm[, c('id', 'Topology')], by = 'id')
        ## signalp
        loc <- which(subloca$X. == 'N') 
        subloca$Signalp[loc] <- 0
        loc <- which(subloca$X. != 'N')
        subloca$Signalp[loc] <- 1
        table(subloca$Signalp)
        ## tmhmm
        loc <- which(subloca$Topology == 'N') 
        subloca$tmhmm[loc] <- 0
        loc <- which(subloca$Topology != 'N')
        subloca$tmhmm[loc] <- 1
        table(subloca$tmhmm)
      }
      
      # --- combine to obtain graph data (all info) & dgi data--- 
      if (anno.info.1 != FALSE) {
        graph_data <- merge(df.reduced.info, subloca, by.y = 'id')
        dgi_data <- graph_data[,c((2:(length(df.No0)+1)), (length(graph_data)-1):length(graph_data))]
        rownames(dgi_data) <- graph_data$id
      } else {
        graph_data <- df.reduced.info
        dgi_data <- graph_data[,c((1:(length(df.reduced.info)-2)))]
        rownames(dgi_data) <- graph_data$id
      }
      
      # --- save data ---
      graphdata_path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
      dgidata_path <- paste("./data/", dataset, "/inputdata/dgidata_", group, '.csv', sep = '')
      write.csv(graph_data, graphdata_path, row.names = FALSE)
      write.csv(dgi_data, dgidata_path, row.names = FALSE)
      
      # temp variable
      mylist[group] <- list(graph_data)
    }
  }
  
  return(mylist)
}

# ========== APPENDIX: Beautify ==========
library(RColorBrewer)
library(paletteer) 
library(ggplot2)
# pallete 

pale_8 <- as.vector(paletteer_d('RColorBrewer::Pastel2'))
pale_9 <- as.vector(paletteer_d('ggprism::pastels'))
pale_11 <- as.vector(paletteer_d('khroma::sunset'))
pale_12 <- as.vector(paletteer_d('RColorBrewer::Set3'))
pale_20 <- as.vector(paletteer_d('ggthemes::Tableau_20'))
pale_25 <- c(pale_12, pale_11, '#85A4FE', 'grey')
pale_29 <- c(pale_25, pale_8)

top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 16),
        legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
mytheme1 <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 16),
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        # axis.title = element_blank(), 
        # legend.title = element_blank(),
        legend.position = 'right', 
        # legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))