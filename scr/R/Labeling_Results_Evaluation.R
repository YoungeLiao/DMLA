# === define function ===
basic_var <- function(dataset, group, SubDataset){
  # --- define path ---
  label_path <- paste("./data/", dataset, '/output_data/result/', SubDataset, '/types.txt', sep = '')
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
  kegg_anno <- read.csv(kegg_anno.path)
  Sub_Swiss <- read.csv(Sub_Swiss.path)
  # matched <- merge(graphdata[,-c((length(graphdata)-4):(length(graphdata)-1))], Sub_Swiss[,1:5], by = 'id', all.x = FALSE, all.y = FALSE)
  # matched_swiss <- merge(graphdata[-c((length(graphdata)-4):(length(graphdata)-1))],
  #                       Sub_Swiss[,c('id', 'SwissProt_ID', 'SwissProt_Description')],
  #                       by = 'id', all = TRUE)
  matched_swiss <- merge(graphdata[-c((length(graphdata)-4):(length(graphdata)-1))], 
                         Sub_Swiss,
                         by = 'id')
  matched_kegg <- merge(matched_swiss, 
                        kegg_anno,
                        by = 'id')
  
  # return var list 
  var_list <- list(label_path, output_path, graphdata.path, fig.path, kegg_anno.path, Sub_Swiss.path,
                   graphdata, labels, 
                   kegg_anno, Sub_Swiss,
                   matched_swiss, matched_kegg)
  return(var_list)
}


Clustering <- function(dataset, group, SubDataset){
  # --- define path ---
  label_path <- paste("./data/", dataset, '/output_data/result/', SubDataset, '/types.txt', sep = '')
  output_path <- paste("./data/", dataset, "/inputdata/", sep = '') # input data for modeling in python
  graphdata.path <- paste("./data/", dataset, "/inputdata/graphdata_", group, '.csv', sep = '')
  fig.path <- paste("./data/", dataset, "/output_data/Figures", sep = '') 
  # --- load data ---
  graphdata <- read.csv(graphdata.path, header = TRUE)
  labels <- paste('Cluster',as.vector(apply(read.table(label_path, sep='\t', header = F)+1, 2, as.factor)))
  graphdata$Cluster <- labels
  
  # --- visualization ---
  library(tidyverse)
  library(cowplot)
  library(ggplot2)
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
  
  ggsave(filename = paste('ClusteringLabeled_', group,'.pdf', sep = ''), width = 6, height = 5, units = 'in', path = fig.path)
  return(p)
}


Loc_FunGene_Clus <- function(GeneName, 
                             dataset, group, SubDataset){
  varlist <- basic_var(dataset, group, SubDataset)
  matched_kegg <- varlist[[length(varlist)]]
  GeneInfoList <- list()
  GeneClusterList <- list()
  for (i in 1:length(GeneName)){
    GeneInfo <- matched_kegg[grepl(GeneName[i], matched_kegg$level3_pathway_name), ] # cluster 3, 5
    GeneCluster <- GeneInfo[c('id', 'Cluster', 'level3_pathway_name')]
    GeneInfoList[GeneName[i]] <- list(GeneInfo)
    GeneClusterList[GeneName[i]] <- list(GeneCluster)
    GeneClusterList$Phototransduction
  }
  mylist <- list()
  mylist['GeneInfoList'] <- GeneInfoList
  mylist['GeneClusterList'] <- GeneClusterList
  return(mylist)
}


Annotation <- function(dataset, group, SubDataset, 
                       GeneSet_bykeg=0,
                       GeneSet_byswiss=0, 
                       GeneSet_bycluster=0){
  "
  return: only save data in output path
  "
  # load data
  var_list <- basic_var(dataset, group, SubDataset)
  matched_kegg <- varlist[[length(varlist)]]
  
  # extract and save gene set
  if (GeneSet_bycluster[[1]] != 0){
    # --- extract functional gene set by cluster ---
    # -- cluster 1 --
    for (i in 1:length(GeneSet_bycluster[[1]])){
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
      print(paste('The total gene number of', GeneSet_byswiss[[i]], 'is:', length(unique(result$id))) ) # 26 genes in total 
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



# === API & conduction DEMO ===
library(dplyr)
dataset <- 'LY_9samples_metatranscriptomics'
group <- 'Yellow_vs_Dark'
SubDataset <- 'YvsD_SubCell_dgi'
# --- Clustering --
p <- Clustering(dataset, group, SubDataset)
p

# --- Loc_FunGene_Clus ---
# locate the cluster that targeted function genes/pathways assigned to 
GeneName <- c('Phototransduction', "Peroxisome", 'Drug metabolism - cytochrome P450', 'Metabolism of xenobiotics by cytochrome P450', 'Longevity regulating pathway - multiple species')
mylist <- Loc_FunGene_Clus(GeneName,
                           dataset, group, SubDataset)
mylist$GeneClusterList

# --- obtain annotation file ---
GeneSet_bycluster <- list('Cluster 1', 'Cluster 3', 'Cluster 4', 'Cluster 6')
GeneSet_byswiss <- list(c('nitrate', 'nitrite'))
GeneSet_bykegg <- list(c('Phototransduction'))
Annotation(dataset, group, SubDataset, GeneSet_bykegg,
           GeneSet_byswiss=0, 
           GeneSet_bycluster=0)

# === visualization backup ===
# --- style 1 ---
ggplot(graphdata, aes(tSNE1,tSNE2, colour = Cluster)) +  # ,colour = degs
  geom_point(alpha = 0.5, size = 2) + 
  # geom_point(colour="#85A4FE", size = 1) + 
  labs(title = '') + 
  scale_colour_manual(values = pale_25) +
  mytheme

# --- style 2 ---
library(tidyverse)
library(cowplot)
library(ggplot2)
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
p
path <- fig.path
ggsave('Yellow_vs_Dark.pdf', p, width = 6, height = 5, units = 'in', path)


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
mytheme1 <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 16),
        legend.position = 'right', 
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))