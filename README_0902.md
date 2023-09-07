# README

## Usage

- Preparation: 
  - group.csv: put in the inputdata directory, with "Samples", "Group" as the column name
- Step1: 
  - Preprocessing.Rmd
  - Refer to the demo (DEMO1, DEMO2, DEMO3), similarly, change the dataset name, and then conduct DESeq analysis, tsne calculation to obtain filtered differentially expressed genes (DEGs) and graph data list. The datasets in the graph data list are also saved in local folder. 
  - Finally, it returned to a graph data list. Graph data list contains multiple datasets based on the provided datasets. Each has a control group and an experimental group. 
- Step2:
  - GenePanelApp.Rmd
  - 

# Draft

- Fig. 2: Preprocessing: 确认是否filter、保存文件如何，
  - OPTIONAL: 是否能导入函数，并在Rmarkdown中展示
  - 结果：filter在umap、tsne中，删除umap仅保留tsne
  - **graphdata: [id + tsne] + normed expression + [degs + cluster label] + subcellular location**
  - tsne中Condition？ 
- Fig.3: DGI labeling results & Evaluation
  - labeling resutls & 之前的函数 之间的差别，并初步复现
    - Clustering: obtain gene panels (labeling and visualization)
    - Loc_FunGene_Clus: locate the cluster that the interested pathways (based on KEGG) and genes (based on Swiss-Prot) that assigned to
      - <u>Swiss-Prot</u>
  - Annotation: Extract the genes and corresponding annotation
    - 问题：Nitrogen metabolism少一个 —— 绘图不影响，暂时跳过
    - **Cluster enrichment**：画出来了，但cluster对应有问题
- Fig. 4: 
  - Annotation: Extract cluster and correponding 
  - GP_Loc_Vis: visualize the location of interested genes in gene panel
  - EnrichmentAna: 

- Fig. 5
  - Obtain_TopoData: obtain data for topological analysis
  - Landmark genes
  - Ref: 
    - GeneExtract: 获取read.csv的file.path —— Cluster+Annotation （无需enrich结果）
    - GenePanelAnalysis: 获取landmark genes
    - GeneNet: 获取基因调控网络的data
  - Additional File 2: ko was utilized to represent each gene





### Directory 

- Inputdata: data that have been processed for other analysis
- output_data: 

### reference: 

- DGI_CoExp.Rproj
- OSProtein.Rproj
  - /Users/yangliao/Documents/GitHub/OSProtein/02_DatasetAnalysis/scr/DEGs_Function.R



### Done

- Preprocessing: data pretreatment. 
  - Output: DESeq results, Graph data info, DGI raw data for graph learning
  - 



一带一路的清华力量：青年融入一带一路人文交流

- 