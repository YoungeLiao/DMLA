## === r VisualizationSet ===

library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(paletteer)

# theme
top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme <- theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 16),
        legend.position = 'none', 
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))

mytheme1 <- theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.text.x = element_text(size = 10, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
        # legend.title = element_blank(),
        legend.title = element_text(size = 16),
        legend.position = 'right', 
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))

# pallete 
Pale_51 <- as.vector(paletteer_d('ggsci::default_igv'))
pale_9_1 <- as.vector(paletteer_d('ggprism::pastels'))
pale_20_2 <- as.vector(paletteer_d('ggthemes::Tableau_20'))
pale_11_2 <- as.vector(paletteer_d('khroma::sunset'))
pale_12_2 <- as.vector(paletteer_d('RColorBrewer::Set3'))
pale_8_1 <- as.vector(paletteer_d('RColorBrewer::Pastel2'))
pale_29 <- c(pale_20_2, pale_9_1)
pale_28 <- c(pale_8_1, pale_11_2, pale_9_1)
pale_31 <- c(pale_8_1, pale_11_2, pale_12_2)
pale_32 <- c(pale_12_2, pale_20_2)
pale_25 <- c(pale_12_2, pale_11_2, '#85A4FE', 'grey')
pale_29 <- c(pale_25, pale_8_1)

pale_10 <- as.vector(paletteer_d('ggsci::default_jco'))
pale_51 <- as.vector(paletteer_d('ggsci::default_igv'))
cluster_pale <- c('#8DD3C7FF','#BEBADAFF',"#FB8072FF", "#FDB462FF")


# === APP1 functions ===
# basic_var <- function(dataset, group, SubDataset){
#   # --- define path ---
#   label_path <- paste("./data/", dataset, '/output_data/', SubDataset, '/generated_data/types.txt', sep = '')
#   output_path <- paste("./data/", dataset, "/inputdata/", sep = '') 
#   graphdata.path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
#   fig.path <- paste("./data/", dataset, "/output_data/Figures", sep = '') 
#   graphdata.labeled.path <- paste("./data/", dataset, "/inputdata/graphdata_", group, 'labeled.csv', sep = '')
#   # annotation
#   kegg_anno.path <- paste("./data/", dataset, "/rawdata/kegg_annotation_all.csv", sep = '') 
#   Sub_Swiss.path <- paste("./data/", dataset, "/inputdata/SubSwissExpre.csv", sep = '') 
#   # --- load data ---
#   graphdata <- read.csv(graphdata.path, header = TRUE)
#   labels <- paste('Cluster',as.vector(apply(read.table(label_path, sep='\t', header = F)+1, 2, as.factor)))
#   graphdata$Cluster <- labels
#   
#   # save labeled data
#   if (file.exists(graphdata.labeled.path)){
#     print(paste("File", graphdata.labeled.path, "exists."))
#   } else {
#     write.csv(graphdata, graphdata.labeled.path, row.names = FALSE)
#   }
#   
#   # load kegg annotation
#   if (file.exists(kegg_anno.path)){
#     kegg_anno <- read.csv(kegg_anno.path)
#     matched_kegg <- merge(graphdata[,-(length(graphdata)-2)], kegg_anno, by = 'id')
#     Sub_Swiss = NA
#     matched_swiss = NA
#   } else {
#     print('Skipped: lack of kegg annotation data.')
#   }
#   
#   if (file.exists(Sub_Swiss.path)){
#     Sub_Swiss <- read.csv(Sub_Swiss.path)
#     matched_swiss <- merge(graphdata[-c((length(graphdata)-4):(length(graphdata)-1))], Sub_Swiss, by = 'id')
#     matched_kegg <- merge(matched_swiss, kegg_anno, by = 'id')
#   }
#   
#   # return var list 
#   var_list <- list(label_path, output_path, graphdata.path, fig.path, kegg_anno.path, 
#                    Sub_Swiss.path,
#                    graphdata, labels, 
#                    kegg_anno, Sub_Swiss,
#                    matched_swiss, matched_kegg)
#   return(var_list)
# }

basic_var <- function(dataset, group, SubDataset){
  # --- define path ---
  label_path <- paste("./data/", dataset, '/output_data/', SubDataset, '/generated_data/types.txt', sep = '')
  output_path <- paste("./data/", dataset, "/inputdata/", sep = '') 
  graphdata.path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
  fig.path <- paste("./data/", dataset, "/output_data/Figures", sep = '') 
  graphdata.labeled.path <- paste("./data/", dataset, "/inputdata/graphdata_", group, 'labeled.csv', sep = '')
  # annotation
  kegg_anno.path <- paste("./data/", dataset, "/rawdata/kegg_annotation_all.csv", sep = '') 
  Sub_Swiss.path <- paste("./data/", dataset, "/inputdata/SubSwissExpre.csv", sep = '') 
  # --- load data ---
  graphdata <- read.csv(graphdata.path, header = TRUE)
  labels <- paste('Cluster',as.vector(apply(read.table(label_path, sep='\t', header = F)+1, 2, as.factor)))
  graphdata$Cluster <- labels
  
  # save labeled data
  if (file.exists(graphdata.labeled.path)){
    print(paste("File", graphdata.labeled.path, "exists."))
  } else {
    write.csv(graphdata, graphdata.labeled.path, row.names = FALSE)
  }
  
  # load kegg annotation
  if (file.exists(kegg_anno.path)){
    kegg_anno <- read.csv(kegg_anno.path)
    matched_kegg <- merge(graphdata, kegg_anno, by = 'id')
    Sub_Swiss = NA
    matched_swiss = NA
  } else {
    print('Skipped: lack of kegg annotation data.')
  }
  
  if (file.exists(Sub_Swiss.path)){
    Sub_Swiss <- read.csv(Sub_Swiss.path)
    matched_swiss <- merge(graphdata[-c((length(graphdata)-4):(length(graphdata)-1))], Sub_Swiss, by = 'id')
    # matched_kegg <- merge(matched_swiss, kegg_anno, by = 'id')
  }
  
  # return var list 
  var_list <- list(label_path, output_path, graphdata.path, fig.path, kegg_anno.path, 
                   Sub_Swiss.path,
                   graphdata, labels, 
                   kegg_anno, Sub_Swiss,
                   matched_swiss, matched_kegg)
  return(var_list)
}


Clustering <- function(dataset, group, SubDataset,
                       reduction = 'tsne'){
  "
  :para reduction: 'tsne' or 'umap'
  "
  var_list <- basic_var(dataset, group, SubDataset)
  # # --- define path ---
  # label_path <- var_list[[1]]
  # output_path <- var_list[[2]]
  # graphdata.path <- var_list[[3]]
  # fig.path <- var_list[[4]]
  
  # --- load data ---
  graphdata <- var_list[[7]]
  labels <- var_list[[8]]
  
  # --- visualization ---
  library(tidyverse)
  library(cowplot)
  library(ggplot2)
  if (reduction == 'tsne'){
    p <- ggplot(graphdata, aes(tSNE1,tSNE2, colour = Cluster)) +
      geom_point(size=2, alpha = 0.8) +
      # stat_ellipse(mapping=aes(group = Cluster,colour = Cluster), #分组边界圈
      #              geom = "path", # other: polygon
      #              linetype = 2,
      #              size=0.6,
      #              alpha=0.5) +
      xlab(NULL) +
      ylab(NULL) +
      theme_cowplot() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        panel.background = element_blank()
      ) +
      scale_colour_manual(values = pale_25)
  }else{
    p <- ggplot(graphdata, aes(umap1, umap2, colour = Cluster)) +
      geom_point(size=2, alpha = 0.8) +
      # stat_ellipse(mapping=aes(group = Cluster,colour = Cluster), #分组边界圈
      #              geom = "path", # other: polygon
      #              linetype = 2,
      #              size=0.6,
      #              alpha=0.5) +
      xlab(NULL) +
      ylab(NULL) +
      theme_cowplot() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        panel.background = element_blank()
      ) +
      scale_colour_manual(values = pale_25)
  }
  
  
  fig.path <- paste("./data/", dataset, "/output_data/Figures", sep = '') 
  ggsave(filename = paste('ClusteringLabeled_', group,'.pdf', sep = ''), width = 6, height = 5, units = 'in', path = fig.path)
  return(p)
}


# Annotate genes 
Annotation <- function(dataset, group, SubDataset, 
                       GeneSet_bykegg=0,
                       GeneSet_byswiss=0, 
                       GeneSet_bycluster=0){
  "
  return: only save data in output path
  "
  # load data
  var_list <- basic_var(dataset, group, SubDataset)
  matched_kegg <- var_list[[length(var_list)]]
  
  # extract and save gene set
  if (GeneSet_bycluster[[1]] != 0){
    # --- extract functional gene set by cluster ---
    # -- cluster 1 --
    for (i in 1:length(GeneSet_bycluster)){
      Cluster <- filter(matched_kegg, matched_kegg$Cluster == GeneSet_bycluster[[i]]) 
      ### save data
      print(paste('The total gene number of', GeneSet_bycluster[[i]], 'is:', length(unique(Cluster$id))) ) # 26 genes in total 
      output_path <- paste("./data/", dataset, "/output_data/R/", SubDataset, '_', GeneSet_bycluster[[i]], '.csv', sep = '') 
      write.csv(Cluster, output_path, row.names = FALSE)
    }
  } else if (GeneSet_byswiss[[1]] != 0){
    # --- extract functional gene set by gene name ---
    for (i in 1:length(GeneSet_byswiss)){
      for (j in 1:length(GeneSet_byswiss[[i]])){
        result.temp <- matched_kegg[grepl(GeneSet_byswiss[[i]][j], matched_kegg$SwissProt_Description), ] 
        if (j == 1){
          result <- result.temp
        } else{
          result <- rbind(result, result.temp)
        }
      }
      print(paste('The total gene number of', GeneSet_byswiss[[i]], 'is:', length(unique(result$id))) ) 
      output_path <- paste("./data/", dataset, "/output_data/R/", SubDataset, '_', GeneSet_byswiss[[i]][j], '.csv', sep = '') 
      write.csv(result, output_path, row.names = FALSE)
    }
  } else if (GeneSet_bykegg[[1]] != 0){
    # --- extract functional gene set by gene name ---
    for (i in 1:length(GeneSet_bykegg)){
      for (j in 1:length(GeneSet_bykegg[[i]])){
        result.temp <- matched_kegg[grepl(GeneSet_bykegg[[i]][j], matched_kegg$level3_pathway_name), ] 
        if (j == 1){
          result <- result.temp
        } else{
          result <- rbind(result, result.temp)
        }
      }
      print(paste('The total gene number of', GeneSet_bykegg[[i]], 'is:', length(unique(result$id))) ) # 26 genes in total 
      output_path <- paste("./data/", dataset, "/output_data/R/", SubDataset, '_', GeneSet_bykegg[[i]][j], '.csv', sep = '') 
      write.csv(result, output_path, row.names = FALSE)
    }
  }
}


# --- define path ---
obtain_pathdata <- function(dataset, group, SubDataset, 
                            Cluster){
  rpkm.path <- paste('./data/', dataset, "/rawdata/RPKM.txt", sep = '')
  Cluster.path <- paste('./data/', dataset, "/output_data/R/", SubDataset, '_', Cluster, '.csv', sep = '')
  group.path <- paste('./data/', dataset, "/inputdata/group.csv", sep = '')
  kegg_anno.path <- paste('./data/', dataset, "/rawdata/kegg_annotation_all.csv", sep = '') 
  Enrichment.path <- paste('./data/', dataset, '/output_data/R/', SubDataset, '_', Cluster, '_Enrichment.csv', sep = '')
  
  # --- load data ---
  Cluster.data <- read.csv(Cluster.path, header = T)
  rpkm <- read.table(rpkm.path, sep='\t', header = T)
  group <- read.csv(group.path)
  
  # --- combine data ---
  rawdata <- merge(rpkm, Cluster.data[,-c(2:10)], by = 'id') # -c(2:length(rpkm))
  samples <- colnames(rawdata)[2:length(rpkm)]
  # print('The samples and groups are as follows:')
  # print(samples)
  # print(group)
  return (rawdata)
}


obtain_pvalue <- function(dataset, ref_group, group_sample.data, pathway, annotation.level){
  
  # prepare data
  group <- group_sample.data
  temp <- data.frame(t(pathway))
  colnames(temp) <- temp[1,]
  name <- colnames(temp)
  temp <- data.frame(temp[-1,])
  
  library(ggpubr)
  temp$group <- group$Group
  temp <- data.frame(apply(temp[,1:(length(temp))], 2, as.numeric))
  temp$group <- group$Group
  
  # initial comparison
  compare_data <- temp[, c(1, length(temp))]
  colnames(compare_data) <- c('Pathway', 'group')
  result <- compare_means(Pathway~group, data = compare_data, ref.group = ref_group, 
                          method = 't.test', 
                          paired = FALSE)
  p_value <- data.frame(t(result$p))
  expre_group.temp <- unique(group$Group)
  loc <- match(ref_group, expre_group.temp)
  expre_group <- expre_group.temp[-loc]
  colnames(p_value) <- expre_group
  
  # delete columns that contain all 0
  feature <- temp[, 1:(length(temp)-1)]
  del <- c() # define vector to store values
  for (i in seq(1, ncol(feature))){ # 
    if(sum(feature[, i])==0){
      print(i)
      del <- append(del, -i)
    }
  }
  if (length(del)){
    feature <- feature[, del]
    temp <- cbind(feature, temp$group)
  }
  
  # loop
  if ((length(temp)-1) != 1) {
    for (i in 2:(length(temp)-1)){ # (length(temp)-1)
      compare_data <- temp[, c(i, length(temp))]
      colnames(compare_data) <- c('Pathway', 'group')
      result <- compare_means(Pathway ~ group, data = compare_data, ref.group = ref_group, 
                              method = 't.test', 
                              paired = FALSE)
      p <- data.frame(t(result$p))
      colnames(p) <- expre_group
      p_value <- rbind(p_value, p)
    }
  }
  
  if (length(del)){
    rownames(p_value) <- name[del]
  }else{
    rownames(p_value) <- name
  }
  
  p_value$pathway <- rownames(p_value)
  
  colnames(p_value) <- c(paste('pvalue', colnames(p_value)[1:(length(p_value)-1)], sep = '_'), annotation.level)
  
  df.pathlevel_p <- merge(pathway, p_value, by = annotation.level)
  
  return(df.pathlevel_p)
}


ObtainStatistic <- function(pathway, df.pathlevel_p, group_sample.data){
  # mean <- apply(pathway[,2:(length(pathway)-1)], 1, mean)
  mean <- apply(df.pathlevel_p[,2:(length(df.pathlevel_p)-2)], 1, mean)
  mean_list <- list()
  sd_list <- list()
  for (i in 1:((length(df.pathlevel_p)-2)/3)){
    sd.temp <- as.data.frame(apply(df.pathlevel_p[,(i*3-1):(i*3+1)], 1, sd))
    mean.temp <- as.data.frame(apply(df.pathlevel_p[,(i*3-1):(i*3+1)], 1, mean))
    sd_list[i] <- sd.temp
    mean_list[i] <- mean.temp
  }
  
  for (i in 2:((length(df.pathlevel_p)-2)/3)){
    df.pathlevel_p[SubDataset] <- mean_list[[i]]/mean_list[[1]]
    df.pathlevel_p[unique(group_sample.data$Group)[i]] <- mean_list[[i]]
  }
  return(df.pathlevel_p)
}


last_chara <- function(string, chara=';'){
  result <- sub(paste0(".*", chara), "", string)
  return (result)
}


enrich <- function(df.pathlevel_p, SubDataset, dataset, Cluster,
                   p_thre, FC_up, FC_down, Expre,
                   group_sample.data, 
                   annotation.level='level3_pathway_name'){
  Path.filtered <- filter(df.pathlevel_p, (df.pathlevel_p[SubDataset] > FC_up | df.pathlevel_p[SubDataset] < FC_down) & 
                            df.pathlevel_p[colnames(df.pathlevel_p[length(df.pathlevel_p)-2])] < p_thre & 
                            df.pathlevel_p[colnames(df.pathlevel_p[length(df.pathlevel_p)])] > Expre) 
  
  Enrichment.path <- paste('./data/', dataset, '/output_data/R/', SubDataset, '_', Cluster, '_Enrichment.csv', sep = '')
  if (file.exists(Enrichment.path)){
    print(paste('File', Enrichment.path, 'Exist.'))
  } else {
    write.csv(Path.filtered, Enrichment.path, row.names = FALSE)
  }
  plot_data <- Path.filtered
  kegg_anno <- read.csv(paste("./data/", dataset, "/rawdata/kegg_annotation_all.csv", sep = '')) 
  
  
  if (annotation.level == 'level3_pathway_name'){
    # match to obtain kegglevel2 labels
    plot_data$KEGGBrite2 <- kegg_anno[match(plot_data$level3_pathway_name, kegg_anno$level3_pathway_name),]$level2_pathway_name
    
    # extract the last level (after ';')
    plot_data$kegglevel3 <- apply(plot_data['level3_pathway_name'], 2, last_chara)
    plot_data$kegglevel2 <- apply(plot_data['KEGGBrite2'], 2, last_chara)
    
    # obtain data for visualization
    start_index <- length(group_sample.data$Samples)+2
    data <- plot_data[, c(1, start_index:(start_index+5))]
    data$Expression <- apply(plot_data[,(start_index-3):(start_index-1)], 1, mean)
    colnames(data) <- c('level3_pathway_name', 'pvalue', 'FoldChange', 'NHQ_mean', 'KEGGBrite2', 'kegglevel3', 'kegglevel2', 'Expression')
    
    # visualization
    library(ggrepel)
    p <- ggplot(data, aes(x = pvalue, y = FoldChange, size = Expression, color = kegglevel2)) +
      geom_point(alpha=0.7)
    p <- p + geom_text_repel(data = data, 
                             aes(pvalue, y = FoldChange, label = kegglevel3),
                             size=4,color="grey20",
                             point.padding = 0.5,hjust = 0.1,
                             segment.color="grey20",
                             segment.size=0.6,
                             segment.alpha=0.8,
                             nudge_x=0,
                             nudge_y=0, 
                             max.overlaps = 10
    ) + 
      labs(x = 'P-value', y = 'Fold change', title = fig.name)
  } else{
    data <- plot_data[,c(1, ((length(plot_data)-2)):length(plot_data))][1:15,]
    colnames(data) <- c('Definition', 'pvalue', 'FoldChange', 'Expression')
    
    # visualization
    library(ggrepel)
    p <- ggplot(data, aes(x = pvalue, y = FoldChange, size = Expression, color = Definition)) +
      geom_point(alpha=0.7)
    p <- p + geom_text_repel(data = data, 
                             aes(pvalue, y = FoldChange, label = Definition),
                             size=4,color="grey20",
                             point.padding = 0.5,hjust = 0.1,
                             segment.color="grey20",
                             segment.size=0.6,
                             segment.alpha=0.8,
                             nudge_x=0,
                             nudge_y=0, 
                             max.overlaps = 10
    ) + 
      labs(x = 'P-value', y = 'Fold change', title = fig.name)
  }

  
  fig <- p + scale_size(range = c(2, 20), name="Expression") +
    scale_color_manual(values = pale_20_2) + #'#006C75', #36CFC8',
    mytheme1
  
  fig.path <- paste("./data/", dataset, "/output_data/Figures", sep = '') 
  ggsave(filename = paste('Enrichment_', SubDataset,'.pdf', sep = ''), width = 10, height = 7.5, units = 'in', path = fig.path)
  
  return(fig)
}


GPsEnrichCompare <- function(dataset, group, SubDataset, k=7){
  '
  :para k: the number of cluster
  '
  # labeling gene panels
  plotdata_cluster.list <- list()
  for (i in 1:k){
    Cluster <- paste('Cluster', i)
    cluster_path <- paste("./data/", dataset, "/output_data/R/", SubDataset, '_', Cluster, '_Enrichment.csv', sep = '') 
    plotdata_cluster.list[[i]] <- read.csv(cluster_path)
    plotdata_cluster.list[[i]]$Cluster <- Cluster
    if (i == 1){
      plotdata_cluster <- as.data.frame(plotdata_cluster.list[[1]])
    } else {
      plotdata_cluster <- rbind(plotdata_cluster, plotdata_cluster.list[[i]])
    }
  }
  
  # prepare plot data
  plotdata_cluster$Expression <- apply(plotdata_cluster[,7:9], 1, mean)
  plotdata <- plotdata_cluster[,-c(2:7)]
  colnames(plotdata) <- c('level3_pathway_name', 'pvalue', 'FoldChange', group, 'Cluster', 'Expression')
  
  # plot
  library(ggrepel)
  p <- ggplot(plotdata,  
              aes(x = pvalue, y = FoldChange, size = Expression, color = Cluster)) +
    geom_point(alpha=0.7)
  p <- p + geom_text_repel(data = plotdata, 
                           aes(pvalue, y = FoldChange,label = level3_pathway_name),
                           size=4,color="grey20",
                           point.padding = 0.5,hjust = 0.1,
                           segment.color="grey20",
                           segment.size=0.6,
                           segment.alpha=0.8,
                           nudge_x=0,
                           nudge_y=0, 
                           max.overlaps = 11
  ) + 
    labs(x = 'P-value', y = 'Fold change', title = group)
  
  p <- p + scale_size(range = c(1, 15), name="Expression") +
    scale_color_manual(values = pale_10) +
    mytheme1
  
  # save fig
  fig.path <- paste("./data/", dataset, "/output_data/Figures", sep = '') 
  ggsave(filename = paste('GPsErichComparison_', group,'.pdf', sep = ''), width = 9, height = 8, units = 'in', path = fig.path)
  
  return(p)
}


Search_GePa <- function(var_list, SearchGene=TRUE, SearchPathway=FALSE){
  matched_kegg <- var_list[[length(var_list)]]
  
  if (SearchGene == TRUE){
    annotation <- 'ko_des'
    # grep gene info list
    kodes <- matched_kegg[grepl(interested.gene, matched_kegg[annotation][,1]), ] 
    
    # extract cluster 
    gene.clu.distri <- table(kodes$Cluster)
    print(paste('Most of', interested.gene, 'related gene belong to:'))
    print(gene.clu.distri[grepl(max(gene.clu.distri), gene.clu.distri)])
    print(paste('The cluster distribution of', interested.gene, 'is:'))
    print(gene.clu.distri)
    
    # obtain cluster 
    temp <- as.data.frame(gene.clu.distri[grepl(max(gene.clu.distri), gene.clu.distri)])
    Cluster <- rownames(temp)
    # print(paste('Gene', interested.gene, 'belongs to', Cluster))
  } else if (SearchPathway == TRUE){
    annotation <- 'level3_pathway_name'
    # grep gene info list
    level3_path <- matched_kegg[grepl(interested.pathway, matched_kegg[annotation][,1]), ] 
    
    # extract cluster 
    path.clu.distri <- table(level3_path$Cluster)
    print(paste('Most of ', interested.pathway, 'related gene belong to:'))
    print(path.clu.distri[grepl(max(path.clu.distri), path.clu.distri)])
    print(paste('The cluster distribution of', interested.pathway, 'is:'))
    print(path.clu.distri)
    
    # obtain cluster 
    temp <- as.data.frame(path.clu.distri[grepl(max(path.clu.distri), path.clu.distri)])
    Cluster <- rownames(temp)
  }
  
  return(Cluster)
}


ObtainCorNet <- function(dataset, SubDataset, Cluster, group, group_sample.data, 
                         gene_name=FALSE, thre.p = 0.05, thre.r = 0.8, annotation.level='KEGG'){
  # obtain filtered pathway
  Filteredcluster.path <- paste("./data/", dataset, "/output_data/R/", SubDataset, '_', Cluster, '_Enrichment.csv', sep = '') 
  Filteredcluster <- read.csv(Filteredcluster.path)
  
  # obtain pathway list & matching & obtain gene list
  rawdata <- obtain_pathdata(dataset, group, SubDataset, Cluster)
  
  # obtain ko descrption as labels
  kegg_anno.path <- paste("./data/", dataset, "/rawdata/kegg_annotation_all.csv", sep = '') 
  kegg_anno <- read.csv(kegg_anno.path)

  # for loop 
  if (annotation.level == 'KEGG'){
    path_list <- Filteredcluster$level3_pathway_name
    
    # initial
    path <- path_list[1]
    gene_list <- rawdata[grepl(path, rawdata$level3_pathway_name), ] 
    for (path in path_list[2:length(path_list)]){
      list <- rawdata[grepl(path, rawdata$level3_pathway_name), ]
      gene_list <- rbind(gene_list, list)
    }
    gene_list <- merge(kegg_anno[c('id')], gene_list, by = 'id')
    
    Col.name <- c('id', 
                  group_sample.data$Samples, 
                  # "SignalP", "tmhmm" , 
                  # "SwissProt_ID", "SwissProt_Description",
                  # "level1_pathway_id", "level2_pathway_id",
                  "level1_pathway_name", "level2_pathway_name", "level3_pathway_id", "level3_pathway_name",
                  'ko', 'ko_name', 'ko_des', 'level3_pathway_name')
    CorCol.name <- c('id', group_sample.data$Samples, 'ko_name')
    Hub_group_raw <- unique(gene_list[Col.name])
    Hub_group <- unique(Hub_group_raw[CorCol.name])
    merge_data <- unique(Hub_group)
  }else{
    path_list <- Filteredcluster$Definition
    
    # initial
    path <- path_list[1]
    gene_list <- rawdata[grepl(path, rawdata$Definition), ] 
    for (path in path_list[2:length(path_list)]){
      list <- rawdata[grepl(path, rawdata$Definition), ]
      gene_list <- rbind(gene_list, list)
    }
    
    CorCol.name <- c('id', group_sample.data$Samples, 'Gene')
    Hub_group <- unique(gene_list[CorCol.name])
    merge_data <- unique(Hub_group)
  }

  
  if (gene_name == TRUE) {
    # v1: merge gene_id and gene description
    gene_names <- data.frame(paste(merge_data$ko_name, ' (', merge_data$id, ')', sep = ''))
    colnames(gene_names) <- 'label'
    corr_data <- cbind(gene_names, merge_data[,2:(length(merge_data)-1)])
    rownames(corr_data) <- corr_data$label
    data <- t(corr_data[, -1])
  } else {
    # v2: only gene id
    corr_data <- merge_data
    rownames(corr_data) <- corr_data$id
    data <- t(corr_data[, -c(1, length(corr_data))])
  }
  
  # calculate correlation matrix
  library(psych)
  occor = corr.test(data,use="pairwise",method="spearman",adjust="fdr",alpha=0.05)
  occor.r = occor$r # 取相关性矩阵R值
  occor.p = occor$p # 取相关性矩阵p值
  occor.r[occor.p>thre.p | abs(occor.r) < thre.r ] = 0
  occor.path <- paste("./data/", dataset, "/output_data/R/", SubDataset, '_', Cluster, '_p', thre.p , '_r', thre.r, '_CorNet.csv', sep = '')
  if (file.exists(occor.path)) {
    print(paste('File', occor.path, 'exists.'))
  } else {
    write.csv(occor.r, occor.path)
  }
  
  return(occor.r)
}


ObtainGraph <- function(dataset, SubDataset, Cluster, group, group_sample.data, igraph, gene_name=TRUE, thre.p = 0.05, thre.r = 0.8){
  
  # === nodes ===
  E(igraph)$corr <- E(igraph)$weight
  E(igraph)$weight <- abs(E(igraph)$weight)
  
  # nodes number check
  print(paste('The total nodes number is:',length(V(igraph)$name))) #或 vcount(igraph)
  
  #节点度（Degree）
  V(igraph)$degree <- degree(igraph)
  # V(igraph)$degree
  
  # --- Visualization: Degree distribution ---
  if (ExploreDegree == TRUE){
    degree_dist <- degree.distribution(igraph)[-1]
    degree_num <- 1:max(V(igraph)$degree)
    par(mfrow = c(1, 2))
    hist(V(igraph)$degree, xlab = 'Degree', ylab = 'Frequency',
         main = 'Degree distribution')
    plot(degree_num, degree_dist, log = 'xy', 
         xlab = 'Log-degree',
         ylab = 'Log-intensity', 
         main = 'Log-log degree distribution')
  }
  
  # --- Exploratory Analysis: Relationship between node degrees and the average degree of its "neighbors" --- 
  if (ExploreNeibor == TRUE){
    neighbor_degree <- graph.knn(igraph, V(igraph))$knn
    plot(V(igraph)$degree, neighbor_degree, log = 'xy',
         xlab = 'Log degree', ylab = 'Log average neighbor degree')
  }
  
  #加权度（Weighted degree）
  V(igraph)$weight_degree <- strength(igraph)
  # V(igraph)$weight_degree
  #接近中心性（Closeness centrality）
  V(igraph)$closeness_centrality <- closeness(igraph)
  # V(igraph)$closeness_centrality
  #介数中心性（Betweenness centrality）
  V(igraph)$betweenness_centrality <- betweenness(igraph)
  # V(igraph)$betweenness_centrality
  #特征向量中心性（Eigenvector centrality）
  V(igraph)$eigenvector_centrality <- evcent(igraph)$vector
  # V(igraph)$eigenvector_centrality
  
  # --- Exploratory Analysis: The relationship between three features describing node centrality ---
  if (ExploreNeibor == TRUE){
    library(car)
    scatter3d(V(igraph)$closeness_centrality, V(igraph)$betweenness_centrality, V(igraph)$eigenvector_centrality, xlab = 'Closeness centrality', ylab = 'Betweenness centrality', zlab = 'Eigenvector centrality', surface = FALSE)
    # Relationship between node degree and node centrality, e.g. the relationship with feature vector centrality
    plot(V(igraph)$degree, V(igraph)$eigenvector_centrality,
         xlab = 'Degree', ylab = 'Eigenvector centrality')
  }
  
  # obtainn node list
  node_list <- data.frame(node_id = V(igraph)$name,
                          degree = V(igraph)$degree,
                          weight_degree = V(igraph)$weight_degree,
                          closeness_centrality = V(igraph)$closeness_centrality,
                          betweenness_centrality = V(igraph)$betweenness_centrality,
                          eigenvector_centrality = V(igraph)$eigenvector_centrality)
  
  
  # save node
  node.path <- paste("./data/", dataset, "/output_data/R/", SubDataset, '_', Cluster, '_p', thre.p , '_r', thre.r, '_Node.csv', sep = '')
  write.table(node_list, node.path, sep = 't', row.names = FALSE, quote = FALSE)
  
  # === Edges ===
  #权重（Weighted）
  E(igraph)$weight
  #边介数中心性（Edge betweenness centrality）
  E(igraph)$betweenness_centrality <- edge.betweenness(igraph)
  E(igraph)$betweenness_centrality
  
  # obtain edge list
  edge <- data.frame(as_edgelist(igraph)) #igraph 
  edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(igraph)$weight,
    correlation = E(igraph)$corr,
    betweenness_centrality = E(igraph)$betweenness_centrality
  )
  
  # save edge
  edge.path <- paste("./data/", dataset, "/output_data/R/", SubDataset, '_', Cluster, '_p', thre.p , '_r', thre.r, '_Edges.csv', sep = '')
  write.table(edge_list, edge.path, sep = 't', row.names = FALSE, quote = FALSE)
  
  graph_list <- list(node_list, edge_list)
  return(graph_list)
}


obtain_landmark <- function(graph_list, dataset, group, SubDataset, 
                            Cluster, group_sample.data,
                            annotation.level='KEGG'){
  # landmark genes: merge nodes list with expression level 
  node_list <- graph_list[[1]]
  edge_list <- graph_list[[2]]
  node_list$id <- node_list$node_id
  var_list <- basic_var(dataset, group, SubDataset)
  matched_kegg <- var_list[[length(var_list)]]
  
  rawdata <- obtain_pathdata(dataset, group, SubDataset, Cluster)
  if (annotation.level == 'KEGG'){
    rawdata <- merge(matched_kegg[c('id')], rawdata, by = 'id')
    
    CorCol.name <- c('id', group_sample.data$Samples, 'ko_name', 'ko_des')
  }else{
    rawdata <- merge(matched_kegg[c('id')], rawdata, by = 'id')
    CorCol.name <- c('id', group_sample.data$Samples, 'Gene', 'Definition')
  }
 
  landmark_info <- unique(rawdata[CorCol.name])
  landmark <- merge(node_list, landmark_info, by = 'id')

  landmark$Mean <- apply(landmark[group_sample.data$Samples], 1, mean)
  landmark.ordered <- landmark[order(landmark$Mean, decreasing = TRUE),]
  return(landmark.ordered)
}