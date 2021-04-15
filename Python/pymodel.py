# %%
%reload_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

import os
from numpy.lib.npyio import load
from scipy.io import loadmat
# from aux import Struct
import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

from mydataset import dataset

device = torch.device('cude' if torch.cuda.is_available() else 'cpu')
print(f'device: {device}')


# %% Load MAT file
fn = 'data_MUA-ONLY_(08-Jan-2021)_bw(5)_fbands(30)_spec(gammatone).mat'
path = './'
Hdry = dataset( drr_idx = 0, unit_number = 1, fn = fn, path = path)
Hdrr = dataset( drr_idx = 4, unit_number = 1, fn = fn, path = path)

# * NORMALIZE
mu = Hdry[:].mean()
std = Hdry[:].std()
Hdry.v = (Hdry.v - mu)/std
Hdrr.v = (Hdrr.v - mu)/std

print(Hdrr)


# %%
train_window = 20   # binwidth * train_window => duration in msec
inout_seq = []      # in/out sequence
for ii in range(Hdry.n_time-train_window):
    seq_ii = Hdrr[ii:ii+train_window]
    label_ii = Hdry[ii+train_window]
    inout_seq.append((seq_ii, label_ii))

test_size = int(0.2*len(inout_seq))
train_seq = inout_seq[:-test_size]
test_seq = inout_seq[-test_size:]       # ! USE JACK-KNIFE !


# # %% Pack data into DATASETs
# test_size = int(0.2 * len(Hdrr))

# trainset = TensorDataset(Hdrr[:-test_size], Hdry[:-test_size])
# testset = TensorDataset(Hdrr[-test_size:], Hdry[-test_size:])

# plt.plot(torch.stack(testset[:],1))
# plt.title('testset')
# plt.xlabel('Samples')
# plt.ylabel('Amp')

# # Create DATALOADERs
# batch_size = 64
# shuffle = False

# train_loader = DataLoader(trainset, batch_size=batch_size, shuffle=shuffle)
# test_loader  = DataLoader(testset, batch_size=batch_size, shuffle=shuffle)


# %% Create RNN model
class RNN(nn.Module):
    def __init__(self, n_input, n_hidden, n_layers, n_output,
        nonlinearity='relu', batch_first=True, dropout=0) -> None:
        super(RNN, self).__init__()
        self.n_input      = n_input         # dimensions of the input tensor
        self.n_hidden     = n_hidden        # num of hidden parameters
        self.n_layers     = n_layers        # num of recurrent RNN layers
        self.nonlinearity = nonlinearity    # {'tanh', 'relu'}
        self.batch_first  = batch_first     # input and output tensors are provided as (batch, seq, feature)
        self.dropout      = dropout 

        self.rnn = nn.RNN(input_size=n_input, hidden_size=n_hidden, num_layers=n_layers,
            nonlinearity=nonlinearity, batch_first=batch_first, dropout=dropout)

        self.fc = nn.Linear(n_hidden, n_output)

    def forward(self, x):
       # (num_layers * num_directions, batch_size, hidden_size)
        self.h0 = torch.randn(self.n_layers, x.shape[0], self.n_hidden)

        # (batch, seq_len, input_size); batch_first == True
        x = x.view(len(x), 1, -1)       #! batch_size == len(x) == train_window

        x, hn = self.rnn(x, self.h0)
        x = self.fc( x[:,-1,:] )          
        return x[-1]


# Create a model
n_input = 1
n_hidden = 100
n_layers = 1
n_output = 1
nonlinearity = 'relu'
batch_first = True
dropout = 0
model = RNN(n_input = n_input, 
            n_hidden = n_hidden, 
            n_layers = n_layers, 
            n_output = n_output,
            nonlinearity = 'relu', 
            batch_first = True, 
            dropout = 0).to(device)

print(model)
# model(train_seq[0][0])  #! *TEST*



# %% 
# *** Training ***
epochs = 100
lr = 1e-3       # learning rate

loss_function = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=lr)

model.train()
loss = []

for ep in range(epochs):
    for k, (seq, y) in enumerate(train_seq):
        optimizer.zero_grad()

        #! Init the hidden states!!
         
        y_pred = model( seq.to(device) )
        loss_k = loss_function(y_pred, y[None])
        loss_k.backward()
        optimizer.step()


    loss.append(loss_k.item())
    if ep%25 == 1:
        print(f'epoch: {ep:4}, loss: {loss[-1]}')

print('- Finished training!')



# %% 
# *** Testing ***
model.eval()
y_pred = []
y_test = []
for k, (seq, yk) in enumerate(test_seq):
    with torch.no_grad():   # ? DO I NEED no_grad() if I have eval() ? 
        yk_pred = model( seq.to(device) )
        y_pred.append(yk_pred)
        y_test.append(yk)

y_pred = np.array(y_pred)
y_test = np.array(y_test)
plt.plot(np.c_[y_pred, y_test, Hdrr[-test_size:]])

CCdry = np.corrcoef(y_pred, y_test)
CCdrr = np.corrcoef(y_pred, Hdrr[-test_size:])
print(f'CCdry: {CCdry[0,1]:8.3f}')
print(f'CCdrr: {CCdrr[0,1]:8.3f}')


# %%

