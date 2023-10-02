# Here is configuration setting.
dataset_name = 'ZJ_12samples_metatranscriptomics'
group = 'Blue_vs_White'
# SubDataset <- 'YvsD_SubCell_dgi'

config = {
    # === DGI: datasets ===
    'data_path': '../data/' + dataset_name + '/inputdata/dgidata_' + group + '.csv', # to do
    'generated_data_path': '../data/' + dataset_name + '/output_data/' + group +'/generated_data/',
    'args_result_path': '../data/' + dataset_name + '/output_data/' + group +'/generated_data/',
    'embedding_data_path': '../data/' + dataset_name + '/output_data/' + group +'/generated_data/',
    'args_model_path': '../data/' + dataset_name + '/output_data/' + group +'/generated_data/',

    # === data preprocessing ===
    'threshold' : 0.04, # 0.4, 4
    'num_feature': 12,

    # === training ===
    'args_n_clusters': 7, # number of clusters desired
    'args_cluster': 10, # run cluster or not 
    'args_num_epoch': 10000, # epoch
    'batch_size': 1,
    'args_DGI': 1,
    'args_lambda_I': 0.8,
    'args_hidden': 256,
    'args_load': 0,
    'args_PCA': 1,
    'args_draw_map': 1,
    'args_diff_gene': '1'
}
