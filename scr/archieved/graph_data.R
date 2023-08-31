# ========== define function ==========
# --- based on python results: labeled matrix ---
get_umap_df <- function(DATAPATH){
  # Import your data
  df <- read.csv(file = DATAPATH, header = TRUE)
  # Generate a new dataset without the 'label' column
  df_umap <- df[, colnames(df) != 'labels' & colnames(df) != 'p_value' & colnames(df) != 'GeneID']
  # Standardize the data
  df_umap <- data.frame(t(apply(df_umap, 1, function(v) {
    (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  })), stringsAsFactors = FALSE)
  # Replace NA values with 0
  df_umap[is.na(df_umap)] <- 0
  # UMAP analysis
  set.seed(0)
  umap <- umap(df_umap, method = 'naive', n_neighbors = 10)
  # Create a data frame for plotting
  df1 <- data.frame(umap$layout)
  df1$label <- rawdata.filtered[,3] # labeling 
  return(df1)
}



# ========== DEMO 1 ==========
library(umap)
library(tidyverse)
library(cowplot)
library(ggplot2)
# --- load data ---
DATAPATH <- './data/ZJ_EES/output_data/ToRdata/labeled_matrix.csv'
df <- get_umap_df(DATAPATH)
# graph_data <- read.csv('./data/ZJ_33samples/inputdata/graph_data.csv', header = TRUE)



rawdata <- read.csv('./data/ZJ_EES/output_data/ToRdata/labeled_matrix.csv', header = TRUE) # all genes have p-value < 0.05 and RPKM > 1/len(counts_df.iloc[:,0])
rawdata$mean <- apply(rawdata[,4:15], 1, mean)
library(dplyr)
# rawdata.filtered <- filter(rawdata, rawdata$mean > 1 & rawdata$p_value < 0.05)
rawdata.filtered <- filter(rawdata, rawdata$p_value < 0.05)
rawdata.norm <- log(rawdata.filtered[,4:15]+1)
# Standardize the data
df_umap <- data.frame(t(apply(df1[,1:2], 1, function(v) {
  (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
})), stringsAsFactors = FALSE)


ggplot(df1, aes(x = X1)) + 
  geom_density(color = 'black', fill = 'gray')



# --- visualization ---

p <- ggplot(graph_data, aes(x = umap1, y = umap2, colour = label)) +
  geom_point(size=0.01) +
  # stat_ellipse(mapping=aes(group = label,colour =label), #分组边界圈
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
  scale_colour_manual(values = c("#0072B5FF","#20854EFF","#8B0000","#E18727FF","#800080","#BC3C29FF"))

p
