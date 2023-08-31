# Here is configuration setting.
dataset_name = 'ZJ_EES'

config = {
    # === DGI: datasets ===
    'data_path': '../data/' + dataset_name + '/inputdata/dgi.csv', # to do
    'generated_data_path': '../data/' + dataset_name + '/output_data/generated_data/',
    'args_result_path': '../data/' + dataset_name + '/output_data/generated_data/',
    'embedding_data_path': '../data/' + dataset_name + '/output_data/generated_data/', # './output/generated_data/',
    'args_model_path': '../data/' + dataset_name + '/output_data/generated_data/', # './output/generated_data/',

    # === data preprocessing ===
    'threshold' : 0.4,
    'num_feature': 18,

    # === training ===
    'args_n_clusters': 6, # number of clusters desired
    'args_cluster': 10, # run cluster or not 
    'args_num_epoch': 20000, # epoch
    'batch_size': 1,
    'args_DGI': 1,
    'args_lambda_I': 0.8,
    'args_hidden': 256,
    'args_load': 0,
    'args_PCA': 1,
    'args_draw_map': 1,
    'args_diff_gene': '1'
}
