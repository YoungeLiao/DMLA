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
    print(paste('The total number of unique DEGs (not unique):',length(unique(DEGs_raw$id))))
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
get_tsne_df <- function(rpkm_path, DEGs.raw, DEGs.filtered.path, 
                        output_path.graphdata, output_path.dgi,
                        anno.sub.path, anno.kegg.path){
  # import data
  rpkm <- read.table(rpkm_path, sep='\t', header = T,row.names = 1)
  rpkm$mean <- apply(rpkm, 1, mean)
  rpkm$id <- rownames(rpkm)
  # DEGs.raw$id <- rownames(DEGs.raw)
  print(paste('The total number of DEGs (not unique):',length(DEGs.raw$id)))
  print(paste('The total number of unique DEGs (not unique):',length(unique(DEGs.raw$id))))
  print(paste('The number of total genes (unique):',length(unique(rpkm$id))))
  
  # merge data
  DEGs <- DEGs.raw[,-c(1,3,4,6)]
  loc <- match(DEGs$id, rpkm$id) # return location in rpkm, so that the abundance rpkm can be matched to DEGs
  rpkm.value <- rpkm[loc,1:(length(rpkm)-1)]
  df.matched <- cbind(DEGs, rpkm.value)
  head(df.matched)
  print(paste('The total number of matched DEGs (not unique):',length(df.matched$id)))
  print(paste('The total number of unique matched DEGs:',length(unique(df.matched$id))))

  # === data for tsne ===
  # head(DEGs.filtered) # checking code
  df <- DEGs.filtered[,-c(1:3)]
  head(df) # checking code
  df.unique <- unique(df[,1:(length(rpkm)-1)])
  # head(df.unique)
  print(paste('The total number of unique DEGs for tsne analysis:', length(df.unique$id)))
  
  if (file.exists(output_path.graphdata)) {
    print(paste("File", output_path.graphdata, "exists."))
    graph_data <- read.csv(output_path.graphdata)
    dgi_data <- read.csv(output_path.dgi)
  } else {
    # Standardize the data
    df_tsne <- data.frame(t(apply(df.unique[,2:(length(rpkm)-1)], 1, function(v) {
      (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    })), stringsAsFactors = FALSE)
    # Replace NA values with 0
    df_tsne[is.na(df_tsne)] <- 0  
    
    # tsne analysis
    set.seed(0)
    library(Rtsne)
    tsne_fpkm = Rtsne(
      unique(df_tsne),
      dims = 2,
      pca = T,
      max_iter = 1000,
      theta = 0.4,
      perplexity = 10,
      verbose = F
    )
    # merge tsne results
    tsne_result = as.data.frame(tsne_fpkm$Y)
    colnames(tsne_result) = c("tSNE1","tSNE2")
    # tsne_result$id <- unique(df.unique[,2:(length(rpkm)-1)])$id # 
    tsne_result$id<- df.unique$id
    # head(tsne_result)
    loc <- match(DEGs.filtered$id, tsne_result$id) # return location in rpkm, so that the abundance rpkm can be matched to DEGs
    tsne_result.values <- tsne_result[loc,1:2]
    # head(tsne_result.values)
    graph_data <- cbind(DEGs.filtered, tsne_result.values)
    # head(graph_data)
    # --- Integrate subcellular information ---
    anno.sub <- read.csv(anno.sub.path)
    graph_data <- cbind(graph_data, anno.sub[,1:3])
    anno.kegg <- read.csv(anno.kegg.path)
    sub <- merge(anno.kegg, anno.sub, by.y = 'id')
    
    # create dataset for dgi 
    dgi_data.tsne <- graph_data[,c('id', 'tSNE1', 'tSNE2')]
    dgi_data.values <- graph_data[,c(5:(length(rpkm)+2))]
    # head(dgi_data.values)
    dgi_data <- unique(cbind(dgi_data.tsne, dgi_data.values))
    
    # --- save data ---
    write.csv(graph_data, output_path.graphdata, row.names = FALSE)
    write.csv(dgi_data, output_path.dgi, row.names = FALSE)
  }
  
  return(graph_data)
}

# ========== DEMO 2 ==========
# umap
library(umap)
library(tidyverse)
library(cowplot)
library(ggplot2)
# tsne
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
output_path.rawDEGs <- paste(output_path, 'DEGs_raw.csv', sep = '')

# --- filtered DEGs ---
rpkm_path <- paste("./data/", dataset, "/rawdata/RPKM.txt", sep = '')
DEGs.filtered.path <- paste("./data/", dataset, "/inputdata/DEGs_filtered.csv", sep = '')

# --- tsne & graph data & dgi---
output_path.dgi <- paste(output_path, 'dgi_tsne.csv', sep = '')
output_path.graphdata <- paste(output_path, 'graphdata_tsne.csv', sep = '')
## annotation
anno.sub.path <- paste("./data/", dataset, "/inputdata/Annotation_Expre.csv", sep = '')
anno.kegg.path <- paste("./data/", dataset, "/inputdata/kegg_annotation_all.csv", sep = '')
## integrate subcellular information 
sub_path <- paste("./data/", dataset, "/inputdata/subloca.csv", sep = '')
# ("./02_DatasetAnalysis/output_data/subloca.csv", row.names = 1)

# obtain DEGs
DEGs.filtered.merged <- obtain_DEGs(rawdata_path, group_path, 
                                    rpkm_path,
                                    output_path.rawDEGs, DEGs.filtered.path)
graph_data <- get_tsne_df(rpkm_path, DEGs.raw, DEGs.filtered.path, output_path.graphdata, output_path.dgi)

# ========= visualization =========
# --- density plot ---
ggplot(graph_data, aes(x = tSNE1)) + 
  geom_density(color = 'black', fill = 'gray')

# ---  raw plot ---
ggplot(graph_data,aes(umap1,umap2)) + 
  geom_point(alpha = 0.4, size = 2) + 
  # geom_point(colour="#85A4FE", size = 1) + 
  labs(title = '') + 
  theme_bw() + 
  scale_colour_manual(values = pale_25) +
  mytheme1

# --- tsne ---
ggplot(graph_data,aes(tSNE1,tSNE2,colour = degs)) + 
  geom_point(alpha = 0.5, size = 1) + 
  # geom_point(colour="#85A4FE", size = 1) + 
  labs(title = '') + 
  theme_bw() + 
  scale_colour_manual(values = pale_25) +
  mytheme1

# ========== APPENDIX ==========

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


