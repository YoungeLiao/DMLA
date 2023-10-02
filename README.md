# DMLA
DLMA: Discover-model-learn-advance cycle for unlocking nature as a code base, i.e. encrypting new biological mechanisms and empowering biotechnology development. 

![WetDryLab](Figures/WetDryLab.png)

This is the repository for Research Article: **From mechanism to application: decrypting light-regulated denitrifying microbiome through geometric deep learning**

## 1 Dev. Environment & Installation

- R  4.2.2
- Python 3.9.17
- Rstudio 2023.06.2+561

For geometric deep learning, please conduct the following command first for installation: 

```
# create virtual environment
conda env create venv_name
conda activate venv_name

# install required packages
pip install -r requirements.txt
```

## 2 Utilization tutorial

### 2.1 Directory 

```
├── data # datasets
│   ├── LY_15samples_metagenomics # Demo dataset 1
│   ├── LY_9samples_metatranscriptomics # Demo dataset 2
│   ├── ZJ_24samples_metagenomics # Demo dataset 3
│   └── ZJ_24samples_metatranscriptomics # Demo dataset 4
├── DMLA.Rproj # R project
├── Figures # Figures for README documentation
│   ├── WetDryLab.png
│   └── Workflow.png
├── GenePanelApp.html # Gene panel apps demenstration in html format
├── GenePanelApp.Rmd # Interactive R markdown demenstration of gene panel apps
├── LICENSE 
├── Preprocessing.Rmd # Data preprocess R markdown demenstrations
├── README_0902.md
├── README.md # Documentation and Instructions
├── requirements.txt # package requirements for installation
└── scr # scripts for geometric deep learning 
    ├── config.py # configuration, i.e. parameter setting 
    ├── main.py # main python file for conducting graph learning 
    ├── Model.py # model
    ├── __pycache__
    ├── R # Supplementary R scripts
    ├── Results_Analysis.py # results analysis
    └── utils.py # additional utilities
```



### 2.2 Data: preparation, processing and output

- data/demo_dataset_name/rawdata/: rawdata, including
  1. **reads.number.txt**: gene expression/abundance represented by reads number
  2. **RPKM.txt**: normalized gene expression/abundance represented by RPKM or FPKM
  3. **kegg_annotation_all.csv**: KEGG annotation
  4. **tmhmm.csv**: subcellular information characterized by whether has a transmembrane domain
  5. **Signal.csv:** subcellular information characterized by whether or not signal peptide

- data/demo_dataset_name/inputdata/: data that have been processed or created manually for further analysis. Files that needed to be provided:
  1. **group.csv**: group information, including 'Samples' and their corresponding 'Group'
  2. **Group_xxx.csv**: xxx represents the experimental group that has been utilized for DESeq, enrichment and other analysis.
- data/demo_dataset_name/outputdata/: output data, including Figures, bioinformatic analysis results and etc. Most of files subjected to this folder are mainly utilized for plotting. 

### 2.3 Usage

- Step 1: Prepare folder and datasets

  - Create your new dataset folder as follow: 

  ```
  ├── inputdata
  ├── output_data
  │   ├── Figures
  │   └── R
  └── rawdata
  ```

  - And then prepare the datasets as described in Section 2.2 

- Step 2: Preprocess data for geometric deep learning 

  - Demonstration: **Preprocessing.Rmd**
    - Refer to the demo (DEMO1, DEMO2, DEMO3), similarly, change the dataset name, and then conduct DESeq analysis, tsne calculation to obtain filtered differentially expressed genes (DEGs) and graph data list. The datasets in the graph data list are also saved in local folder. 
    - Finally, it returned to a graph data list. Graph data list contains multiple datasets based on the provided datasets. Each has a control group and an experimental group. 
  - Output: Detailed explaination please refer to Preprocessing.Rmd. All the processed data are saved in the 'inputdata' directory, including: 
    - DEGs: Differentially expressed genes, including DESeq raw results (DEGs_raw.csv) and filtered DEGs (DEGs_filtered.csv)
    - dgidata: such as 'dgidata_Yellow_vs_Dark.csv'. This kind of data is to used for geometric deep learning. Yellow_vs_Dark represents the experimental and control groups, indicating we are going to explore how yellow light can influence or regulate the biological system.
    - graphdata: information of targeted subdataset (e.g. blue vs dark), including genes, expression level of the control and experimental samples, DEGs label, subcellular raw data.

- Step 3: Conduct geometric deep learning 

  - Demenstration: **scr/main.py**
    - Step 3.1: In script 'config.py', change the 'dataset_name' to the name of your dataset, and 'group' to the targeted comparison subdataset.
    - Step 3.2: run 'main.py'
  - Output: The results will be saved in the 'output_data/subdataset/'. 'types.txt' is the labeling results of gene panels that can be matched to the graphdata or dgidata. 

- Step 4: Gene panel related apps

  - Demonstration: **GenePanelApp.Rmd**
  - There are 4 apps as follows: 
    - APP1: 2D-visualization of gene panels
    - APP2: Gene panel enrichment analysis for unlocking unknown potentials of microbiomes
    - APP3: Gene Panel comparison to identify functional gene panels
    - APP4: Topology Network with landmark genes for developing regulatory strategy
  - Output: Please refer to the demonstration and explainations in GenePanelApp.Rmd.

## 3 Case demonstrations 

See: [Preprocessing.html](https://github.com/YoungeLiao/DMLA/blob/main/Preprocessing.html) for data preprocessing, and [GenePanelApp.html](https://github.com/YoungeLiao/DMLA/blob/main/GenePanelApp.html) for making the most of gene panels applications toolkits to explore your natural code base, i.e. the microbiome, and unlock the natural potential for science discovery and new biotechnology development. 

## 4 Citing DMLA

```
@article{liao2023DMLA,
title={From mechanism to application: decrypting light-regulated denitrifying microbiome through geometric deep learning},
author={Yang Liao, Jiyong Bian, Jing Zhao, Ziwei Zhang, Siqi Xu, Yijian Qin, Shiyu Miao, Rui Li, Ruiping Liu*, Meng Zhang, Wenwu Zhu, Huijuan Liu, Jiuhui Qu.},
journal={iMeta},
year={2023},
publisher={Tsinghua University}
}
```

## 5 Contact

If you have any questions or feedback, feel free to reach Yang Liao: liaoy21@mails.tsinghua.edu.cn 
