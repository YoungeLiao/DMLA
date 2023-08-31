
# data name
dataset <- 'LY_18samples_16S'
label_path <- paste("./data/", dataset, "/output_data/generated_data/types.txt", sep = '')
output_path <- paste("./data/", dataset, "/inputdata/", sep = '')
graphdata <- paste(output_path, 'graphdata_umap.csv', sep = '')

# load data
graphdata <- read.csv(graphdata, header = TRUE)
labels <- apply(read.table(label_path, sep='\t', header = F), 2, as.factor)
labeled_data <- cbind(labels, graphdata)

# --- umap ---
p <- ggplot(labeled_data, aes(x = umap1, y = umap2, colour = V1)) +
  geom_point(size=2, alpha = 0.8) +
  stat_ellipse(mapping=aes(group = V1,colour = V1), #分组边界圈
               geom = "path", # other: polygon
               linetype = 2,
               size=0.6,
               alpha=0.5) +
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