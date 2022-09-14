from __future__ import print_function
import argparse
import torch
import torch.nn.functional as F
from torch import nn, optim
from torchvision import datasets, transforms
from tqdm import tqdm
from scipy.io import loadmat, savemat
from torch.autograd import Variable
import numpy as np
from torch.utils.data import TensorDataset, DataLoader
import torch.utils.data as utils
from torchvision.utils import save_image
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing



chrs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X']

def runEachChrs(i):

    class MyConvTranspose2d(nn.Module):
        def __init__(self, conv, output_size):
            super(MyConvTranspose2d, self).__init__()
            self.output_size = output_size
            self.conv = conv
            
        def forward(self, x):
            x = self.conv(x, output_size=self.output_size)
            return x

    class UnFlatten(nn.Module):
        def forward(self, input, size=64):
            return input.view(input.size(0), size, out_shape, out_shape)

    class VAE(nn.Module):
        def __init__(self):
            super(VAE, self).__init__()
            self.invConv1 = nn.ConvTranspose2d(64, 32, kernel_size=4, stride=2, bias=False)
            self.invConv2 = nn.ConvTranspose2d(32, 1, kernel_size=4, stride=2, bias=False)

            self.encoder = nn.Sequential(
                nn.Conv2d(1, 32, kernel_size=4, stride=2, bias=False),
                nn.ReLU(True),
                nn.Conv2d(32, 64, kernel_size=4, stride=2, bias=False),
                nn.ReLU(True),
                nn.Flatten()
            )

            self.fc1 = nn.Linear(h_dim, 20)
            self.fc2 = nn.Linear(h_dim, 20)
            self.fc3 = nn.Linear(20, h_dim)
            
            self.decoder = nn.Sequential(
                UnFlatten(),
                MyConvTranspose2d(self.invConv1, output_size=(out_shape_internal, out_shape_internal)),
                nn.ReLU(True),
                MyConvTranspose2d(self.invConv2, output_size=(tensor_shape[1], tensor_shape[1])),
                nn.Softplus(),
            )
            
        def reparameterize(self, mu, logvar):
            std = logvar.mul(0.5).exp_()
            # return torch.normal(mu, std)
            esp = torch.randn(*mu.size())
            z = mu + std * esp
            return z
        
        def bottleneck(self, h):
            mu, logvar = self.fc1(h), self.fc2(h)
            z = self.reparameterize(mu, logvar)
            return z, mu, logvar

        def forward(self, x):
            h = self.encoder(x)
            z, mu, logvar = self.bottleneck(h)
            t = self.fc3(z)
            zz = self.decoder(t)
            return zz, mu, logvar, z




    def loss_function(recon_x, x, mu, logvar):
        BCE = F.poisson_nll_loss(recon_x , x, reduction='sum', log_input=False)
        KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
        return BCE + KLD



    def train(epoch):
        model.train()
        train_loss = 0
        for _, (data) in enumerate(train_loader):
            data = data[0]
            data = data.to(device)
            optimizer.zero_grad()
            recon_batch, mu, logvar, _ = model(data)
            loss = loss_function(recon_batch, data, mu, logvar)
            loss.backward()
            train_loss += loss.item()
            optimizer.step()
        print('====> Epoch: {} Average loss: {:.4f}'.format(
            epoch, train_loss / len(train_loader.dataset)))

    chr_i = str(i)
    ###Load data set
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    tensor = np.load('./input_folder/' + chr_i + '.npy')

    torch.autograd.set_detect_anomaly(True) 
    a = []
    for i in range(9230):
        a.append(np.sum(tensor[i,:,:]))
    
    for i in range(9230):
        tensor[i,:,:] = tensor[i,:,:] / np.sum(tensor[i,:,:]) * np.mean(a)
    
    tensor_shape = tensor.shape
    def get_dims(shape, kernel_size, stride):
        return (shape - kernel_size) / stride + 1
    out_shape_internal = int(get_dims(tensor_shape[1], 4, 2))
    out_shape = int(get_dims(out_shape_internal, 4, 2))
    h_dim = 64 * out_shape * out_shape

    net_data = np.expand_dims(tensor, axis=1)
    batch_size = 10
    net_data = torch.Tensor(net_data)
    net_data = TensorDataset(net_data)
    train_loader = utils.DataLoader(net_data, batch_size) 

    model = VAE().to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)

    for epoch in range(20):
        # if epoch < 5:
        #     model.fc11.requires_grad = False
        #     model.fc21.requires_grad = False
        # if epoch >= 5:
        #     model.fc11.requires_grad = True
        #     model.fc21.requires_grad = True
        train(epoch)

    train_loader = utils.DataLoader(net_data, 9230) 
    for batch_idx, (data) in enumerate(train_loader):
        net = data[0]
    recon_batch, mu, logvar, z = model(net)
    np.savetxt('./output_folder/latent_' + chr_i + '.csv', z.detach().numpy(), delimiter=",")

results = Parallel(n_jobs=20)(delayed(runEachChrs)(i) for i in chrs)

