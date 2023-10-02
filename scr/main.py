import numpy as np
import pandas as pd
import sys, os
sys.path.append('./scr')
sys.path.append('..')
import pickle

import scipy
from scipy import sparse
import scipy.linalg

from config import config, dataset_name, group
from utils import getNodes, getEdges, get_graph, getEdges_cuda
import utils
import Model
import Results_Analysis

import scanorama
import openpyxl
from matplotlib import pyplot as plt

# import cupy as cp
import torch
from torch_geometric.data import Data, DataLoader

def set_config(DATA_NAME, Subdataset, config=config, n_clusters = 7, thre = 0.4, num_epoch = 20000):
    dataset_name = DATA_NAME
    group = Subdataset
    config['data_path'] = '../data/' + dataset_name + '/inputdata/dgidata_' + group + '.csv' # to do
    config['generated_data_path'] = '../data/' + dataset_name + '/output_data/' + group +'/generated_data/'
    config['args_result_path'] = '../data/' + dataset_name + '/output_data/' + group +'/generated_data/'
    config['embedding_data_path'] = '../data/' + dataset_name + '/output_data/' + group +'/generated_data/'
    config['args_model_path'] = '../data/' + dataset_name + '/output_data/' + group +'/generated_data/'

    config['args_n_clusters'] = n_clusters
    config['args_num_epoch'] = num_epoch # epoch
    config['threshold'] = thre # threshold for edges 


def main():
    import numpy as np
    import pandas as pd
    import sys, os
    from config import config, dataset_name, group
    from utils import getNodes, getEdges, get_graph, getEdges_cuda
    import utils
    import Model
    import Results_Analysis

    import scanorama
    import openpyxl
    from matplotlib import pyplot as plt

    # import cupy as cp
    import torch
    from torch_geometric.data import Data, DataLoader
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = '0,1'
    device_ids = [0, 1]


    # Get nodes features
    Node_features = getNodes(config)


    # get edges
    if torch.cuda.is_available():
        print('cuda is available')
        getEdges_cuda(config) # Edge data will be generated. 
    else:
        print('cuda is not available')
        getEdges(config)


    # Get graph data
    X_data, adj = utils.get_GraphData(config)
    import torch
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    X_data, adj = utils.get_GraphData(config)


    # Trainning 
    if config['args_DGI']:
        print("-----------Deep Graph Infomax-------------")
        data_list = utils.get_graph(adj, X_data)
        data_loader = DataLoader(data_list, batch_size= config['batch_size'])
        np.random.seed(0)
        DGI_model = Model.train_DGI(config, data_loader=data_loader)
        for data in data_loader:
            data.to(device)
            X_embedding, _, _ = DGI_model[0](data)
            X_embedding = X_embedding.cpu().detach().numpy()
            if not os.path.exists(config['generated_data_path']):
                os.makedirs(config['generated_data_path'])
            X_embedding_filename =  config['generated_data_path'] + 'lambdaI' + str(config['args_lambda_I']) + '_epoch' + str(config['args_num_epoch']) + '_Embed_X.npy'
            np.save(X_embedding_filename, X_embedding)


    # plot loss 
    def plot(DGI_model):
        all_loss_list = DGI_model[2]
        plt.plot(range(len(all_loss_list)), all_loss_list, 'b')  
        plt.savefig(config['args_result_path'] + 'lossplot.pdf')
        df_all_loss_list = pd.DataFrame(all_loss_list)
        df_all_loss_list.to_csv(config['args_result_path'] +'/loss.csv')
        np.savetxt(config['args_result_path'] +'/loss.txt', np.array(all_loss_list), fmt='%.4f', delimiter='\t') 
    # plot(DGI_model)


    # obtain embedding
    for data in data_loader:
        data.to(device)
        X_embedding, _, _ = DGI_model[0](data)
        X_embedding = X_embedding.cpu().detach().numpy()
        X_embedding_filename =  config['generated_data_path'] + 'lambdaI' + str(config['args_lambda_I']) + '_epoch' + str(config['args_num_epoch']) + '_Embed_X.npy'
        np.save(X_embedding_filename, X_embedding)


    # Clustering 
    if config['args_cluster']:
        print("-----------Clustering-------------")
        X_embedding_filename =  config['generated_data_path'] +'lambdaI' + str(config['args_lambda_I']) + '_epoch' + str(config['args_num_epoch']) + '_Embed_X.npy'
        X_embedding = np.load(X_embedding_filename)
        if config['args_PCA']:
            X_embedding = Results_Analysis.PCA_process(X_embedding, nps=30)

    # continue the function above
    if config['args_cluster']:   
        print('Shape of data to cluster:', X_embedding.shape)
        cluster_labels, score = Results_Analysis.Kmeans_cluster(X_embedding, config['args_n_clusters']) 
        config['args_n_clusters'] = cluster_labels.max()+1
        all_data = cluster_labels
        generated_data_path = config['generated_data_path']
        if not os.path.exists(config['args_result_path']):
            os.makedirs(config['args_result_path'])
        np.savetxt(config['args_result_path'] + '/types.txt', np.array(all_data), fmt='%3d', delimiter='\t')


    # save clustering labels and print SCI score (Evaluation)
    np.savetxt(config['args_result_path'] + '/types.txt', np.array(all_data), fmt='%3d', delimiter='\t')
    print('SCI score (Clustering quality) is:', score)

if __name__ == '__main__':
    main()