import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, ChebConv, GATConv, DeepGraphInfomax, global_mean_pool, global_max_pool  # noqa
from torch_geometric.data import Data, DataLoader
import matplotlib.pyplot as plt

import hiddenlayer as hl
import time

class Encoder(nn.Module):
    def __init__(self, in_channels, hidden_channels):
        super(Encoder, self).__init__()
        self.conv = GCNConv(in_channels, hidden_channels)
        self.conv_2 = GCNConv(hidden_channels, hidden_channels)
        self.conv_3 = GCNConv(hidden_channels, hidden_channels)
        self.conv_4 = GCNConv(hidden_channels, hidden_channels)
        
        self.prelu = nn.PReLU(hidden_channels)
        
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_attr
        x = self.conv(x, edge_index, edge_weight=edge_weight)
        x = self.conv_2(x, edge_index, edge_weight=edge_weight)
        x = self.conv_3(x, edge_index, edge_weight=edge_weight)
        x = self.conv_4(x, edge_index, edge_weight=edge_weight)
        x = self.prelu(x)

        return x

class my_data():
    def __init__(self, x, edge_index, edge_attr):
        self.x = x
        self.edge_index = edge_index
        self.edge_attr = edge_attr

def corruption(data):
    x = data.x[torch.randperm(data.x.size(0))]
    return my_data(x, data.edge_index, data.edge_attr)

def train_DGI(config, data_loader):
    import numpy as np
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    DGI_model = DeepGraphInfomax(
        hidden_channels=config['args_hidden'],
        encoder=Encoder(in_channels = config['num_feature'], hidden_channels=config['args_hidden']),
        summary=lambda z, *args, **kwargs: torch.sigmoid(z.mean(dim=0)),
        corruption=corruption).to(device)
    DGI_optimizer = torch.optim.Adam(DGI_model.parameters(), lr=1e-6)
    if config['args_load']:
        DGI_filename = config['generated_data_path']+'DGI_lambdaI_' + str(args_lambda_I) + '_epoch' + str(config['args_num_epoch']) + '.pth.tar'
        DGI_model.load_state_dict(torch.load(DGI_filename))
    else:
        import datetime
        start_time = datetime.datetime.now()
        DGI_loss_list = []
        for epoch in range(config['args_num_epoch']):
            DGI_model.train()
            DGI_optimizer.zero_grad()   
            
            DGI_all_loss = []    
            for step, data in enumerate(data_loader):
                data = data.to(device)
                pos_z, neg_z, summary = DGI_model(data=data)

                DGI_loss = DGI_model.loss(pos_z, neg_z, summary)
                DGI_loss.backward()
                DGI_all_loss.append(DGI_loss.item())
                DGI_optimizer.step()
                

            if ((epoch+1)%100) == 0:
                print('Epoch: {:03d}, Loss: {:.4f}'.format(epoch+1, np.mean(DGI_all_loss)))
                DGI_loss_list.append(np.mean(DGI_all_loss))

        end_time = datetime.datetime.now()
        DGI_filename =  config['generated_data_path']+'DGI_lambdaI_' + str(config['args_lambda_I']) + '_epoch' + str(config['args_num_epoch']) + '.pth.tar'
        torch.save(DGI_model.state_dict(), DGI_filename)
        print('Training time in seconds: ', (end_time-start_time).seconds)
                
    return DGI_model, DGI_all_loss, DGI_loss_list
