# ========== define function ==========
## === DEGs ===
obtain_DEGs <- function(rawdata_path, group_path, 
                        rpkm_path,
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
    DEGs_raw <- read.csv(file_path)
  } else {
    # === Filter Differential expressed genes (DEGs) ===
    # --- threshold 1: padj < 0.05, |log2 fold change| > 1 --- 
    print(paste('The total number of DEGs (not unique):',length(DEGs_raw$id)))
    print(paste('The total number of unique DEGs (not unique):', length(unique(DEGs_raw$id))))
    DEGs.filtered <- subset(DEGs_raw,  padj < 0.05 & abs(log2FoldChange) > 1) ## padj
    print(paste('The total number of deseq filtered DEGs:',length(DEGs.filtered$id))) 
    print(paste('The total number of unique deseq filtered DEGs:', length(unique(DEGs.filtered$id)))) 
    
    # --- threshold 2: mean > 1/length ---
    # obtian expression levels (fpkm)
    rpkm <- read.table(rpkm_path, sep='\t', header = T, row.names = 1)
    rpkm$mean <- apply(rpkm, 1, mean)
    rpkm$id <- rownames(rpkm)
    print(paste('The number of total genes (unique):',length(unique(rpkm$id))))
    loc <- rpkm$mean >= thre
    thre <- 1 # 1/(length(rpkm)-2)
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

## === TSNE ===
DEGs.filtered <- DEGs.filtered.merged
get_tsne_df <- function(DEGs.filtered, dataset,
                        group_path,
                        anno.signal, anno.tmhmm){
  # # import data
  # rpkm <- read.table(rpkm_path, sep='\t', header = T,row.names = 1)
  # rpkm$mean <- apply(rpkm, 1, mean)
  # rpkm$id <- rownames(rpkm)
  # # DEGs.raw$id <- rownames(DEGs.raw)
  # print(paste('The total number of DEGs (not unique):',length(DEGs.raw$id)))
  # print(paste('The total number of unique DEGs (not unique):',length(unique(DEGs.raw$id))))
  # print(paste('The number of total genes (unique):',length(unique(rpkm$id))))
  # 
  # # merge data
  # DEGs <- DEGs.raw[,-c(1,3,4,6)]
  # loc <- match(DEGs$id, rpkm$id) # return location in rpkm, so that the abundance rpkm can be matched to DEGs
  # rpkm.value <- rpkm[loc,1:(length(rpkm)-1)]
  # df.matched <- cbind(DEGs, rpkm.value)
  # head(df.matched)
  # print(paste('The total number of matched DEGs (not unique):',length(df.matched$id)))
  # print(paste('The total number of unique matched DEGs:',length(unique(df.matched$id))))

  # === data for tsne ===
  # head(DEGs.filtered) # checking code
  df <- unique(DEGs.filtered[,-c(2:7, length(DEGs.filtered))])
  head(df) # checking code
  
  # normalize the expression levels 
  # normalization
  df.norm <- log(df[,-(1:2)] + 1, 10)
  df.norm$id <- df$id
  df.norm$degs <- df$degs
  
  # --- grouping to different datasets ---
  DEG_group <- as.vector(as.data.frame(table(DEGs.filtered$degs))$Var1)
  # load conditions information
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
      tsne.temp <- subset(df.norm, df.norm$degs == group)
      tsne.raw <- tsne.temp[,1:(length(tsne.temp)-2)]
      print(paste('The total number of unique DEGs for tsne analysis of', group, 'dataset is:', length(tsne.temp$id)))
      
      condition <- as.vector(
        as.data.frame(table(conditions$Group[-(1:3)]))$Var1[i]
      )
      samples <- subset(conditions, Group == condition | Group == conditions$Group[1])$Samples
      df_tsne <- tsne.raw[,samples]
      
      # --- tsne analysis ---
      library(Rtsne)
      set.seed(0)
      tsne_fpkm = Rtsne(
        unique(df_tsne),
        dims = 2,
        pca = T,
        max_iter = 1000,
        theta = 0.4,
        perplexity = 20,
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
      # load SignalP and tmhmm
      sub.signalp <- read.csv(anno.signal, header = T)
      sub.tmhmm <- read.csv(anno.tmhmm, header = T)
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
      
      # --- combine to obtain graph data (all info) & dgi data--- 
      tsne_info <- cbind(tsne_result, df_tsne, tsne.temp[,c('degs', 'id')])
      graph_data <- merge(tsne_info, subloca, by.y = 'id')
      head(graph_data)
      dgi_data <- graph_data[,c((2:(length(df_tsne)+2)), (length(graph_data)-1):length(graph_data))]
      rownames(dgi_data) <- graph_data$id
      
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

# ========== DEMO ==========
library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(paletteer)

# === import data path ===
# --- DESeq ---
dataset <- 'LY_9samples_metatranscriptomics'
rawdata_path <- paste("./data/", dataset, "/rawdata/reads_number.txt", sep = '')
group_path <- paste("./data/", dataset, "/inputdata/group.csv", sep = '')
output_path <- paste("./data/", dataset, "/inputdata/", sep = '') # input data for modeling in python
# output_path.rawDEGs <- paste(output_path, 'DEGs_raw.csv', sep = '')

# --- filtered DEGs ---
rpkm_path <- paste("./data/", dataset, "/rawdata/RPKM.txt", sep = '')
DEGs.filtered.path <- paste("./data/", dataset, "/inputdata/DEGs_filtered.csv", sep = '')

# --- tsne & graph data & dgi---
output_path.dgi <- paste(output_path, 'dgi_tsne.csv', sep = '')
output_path.graphdata <- paste(output_path, 'graphdata_tsne.csv', sep = '')
## annotation
anno.signal <- paste("./data/", dataset, "/rawdata/Signal.csv", sep = '')
anno.tmhmm <- paste("./data/", dataset, "/rawdata/tmhmm.csv", sep = '')

# === conduct analysis ===
# obtain DEGs
DEGs.filtered <- obtain_DEGs(rawdata_path, group_path, 
                                    rpkm_path,
                                    output_path.rawDEGs, DEGs.filtered.path)
# obtain graph data
graph_data_list <- get_tsne_df(DEGs.filtered, dataset,
                               group_path,
                               anno.signal, anno.tmhmm)


# ========= visualization toolkits =========
# --- density plot ---
ggplot(df.norm, aes(x = Blue1)) + 
  geom_density(color = 'black', fill = 'gray')

# ---  tsne plot ---
graph_data <- graph_data_list$Blue_vs_Dark
ggplot(graph_data,aes(tSNE1,tSNE2, colour = as.factor(tmhmm))) +  # ,colour = degs
  geom_point(alpha = 0.5, size = 1) + 
  # geom_point(colour="#85A4FE", size = 1) + 
  labs(title = '') + 
  theme_bw() + 
  scale_colour_manual(values = pale_25) +
  mytheme1

# ========== APPENDIX: Beautify ==========
library(RColorBrewer)
library(paletteer) 
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


