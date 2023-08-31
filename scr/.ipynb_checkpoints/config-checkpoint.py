# Here is configuration setting.
dataset_name = 'ZJ_EES'

config = {
    # === DESeq: datasets & files ===
    'dataset_name' : dataset_name,
    'DATAPATH' : '../data/' + dataset_name + '/rawdata/reads_number.txt', 
    'FEATUREPATH' : '../data/' + dataset_name + '/inputdata/group.csv',
    'OUTPUT_PATH' : '../data/' + dataset_name + '/output_data/',
    'TEMP_FILENAME' : 'counts_df.csv',
    'RPKM_PATH' : '../data/' + dataset_name + '/rawdata/RPKM.txt',

    # === DESeq: hyperparameters ===
    'condition' : 'Group',
    'ref' : 'Control',
    'factors' : 'Group',

    # === DGI: datasets ===
    'data_path': '../data/' + dataset_name + '/inputdata/graph_data.csv', # './data/graph_data_AllDEGs.csv',
    'generated_data_path': '../data/' + dataset_name + '/output_data/generated_data/', # to do
    'args_result_path': '../data/' + dataset_name + '/output/results/',
    
    # 'args_embedding_data_path': './output/generated_data/', # './output/generated_data/',
    # 'args_data_path': './output/generated_data/', # './output/generated_data/',
    # 'args_model_path': './output/generated_data/', # './output/generated_data/',
    
    # 'ROOT_OUTPUT_DATA_PATH' : '../data/' + dataset_name + './output/generated_data/',
    # 'ROOT_RESULTS_PATH': '../data/' + dataset_name + './output/results/', 
    
    # 'edge_data_path': 'generated_data_path/' + 'edges.npy', # to do

    # === data preprocessing ===
    'threshold' : 0.4,
    'num_feature': 9,

    # === training ===
    'batch_size': 1,

    'args_DGI': 1,
    'args_lambda_I': 0.8,
    'args_n_clusters': 24, # number of clusters desired
    'args_cluster': 10, # run cluster or not 
    'args_num_epoch': 10000,
    'args_hidden': 256,
    'args_load': 0,
    'args_PCA': 1,
    'args_draw_map': 1,
    'args_diff_gene': '1'
}

# config: {'INPUT_SHAPE': (450, 10),
#           'ENCODER_SHAPE': [512, 256],
#           'DECODER_SHAPE': [256, 512],
#           'ACTIVATION': 'relu',
#           'LAST_ACTIVATION': 'softmax',
#           'DROPOUT': 0,
#           'LATENT_SCALE': 5,
#           'OPTIMIZER': 'Adam',
#           'BATCH_SIZE': 64,
#           'EPOCHS': 300,
#           'STEPS_PER_EPOCH': None,
#           'VALIDATION_SPLIT': 0.2,
#           'VALIDATION_STEPS': 10,
#           'LATENT_OFFSET': 10,
#           'DECODER_BIAS': 'last',
#           'DECODER_REGULARIZER': 'var_l1',
#           'DECODER_REGULARIZER_INITIAL': 0.0001,
#           'DECODER_RELU_THRESH': 0,
#           'BASE_LOSS': 'mse',
#           'DECODER_BN': False,
#           'CB_MONITOR': 'val_loss',
#           'CB_LR_USE': True,
#           'CB_LR_FACTOR': 0.2,
#           'CB_LR_PATIENCE': 15,
#           'CB_LR_MIN_DELTA': 1e-8,
#           'CB_ES_USE': True,
#           'CB_ES_PATIENCE': 30,
#           'CB_ES_MIN_DELTA': 0,
#           'MULTI_GPU': False
#           }