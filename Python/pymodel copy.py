# %%
# %reload_ext autoreload
# %autoreload 2
# %reset -f 
from IPython import get_ipython
ipython = get_ipython()

if '__IPYTHON__' in globals():
    ipython.magic('load_ext autoreload')
    ipython.magic('autoreload 2')
    ipython.magic('reset -f')


# %%    
from myaux import whitening
import matplotlib.pyplot as plt
import matplotlib
# from numpy.core import tests
# from sklearn.utils import shuffle
matplotlib.style.use('ggplot')

import os
from numpy.lib.npyio import load
from scipy.io import loadmat
# from aux import Struct
import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

from sklearn.model_selection import train_test_split

from mydataset import dataset

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f'device: {device}')



# %% Load from MAT file
fn = 'data_MUA-ONLY_(08-Jan-2021)_bw(5)_fbands(30)_spec(gammatone).mat'
path = './'

seq_len = 20                 # samples   #* hyper-parameter 
unit_numbers = (0,1)                    #* NUMBER of UNITS for analysis 

dry_data = dataset( 
    seq_len = seq_len,      # binwidth * seq_len => duration in msec    
    drr_idx = 0, 
    unit_numbers = unit_numbers, 
    fn = fn, 
    path = path)

fn_strf = 'STRF_MUA-ONLY_bw(5)_fbands(30)_spec(gammatone).mat'
drr_data = dataset(
    seq_len = seq_len,      # binwidth * seq_len => duration in msec    
    drr_idx = 0, 
    unit_numbers = unit_numbers, 
    fn = fn_strf, 
    path = path)

print(drr_data)


 
# %% Split the data into train/test sets
test_size = 0.1
batch_size = 64
shuffle = False

train_idx, test_idx = train_test_split(range(len(dry_data)), test_size=test_size, 
    shuffle=shuffle )  #, random_state=42)


# T0 = test_idx   # ! ### DEBUG ###

# # ! ### DEBUG ###
# train_idx = np.arange(600, len(dry_data))
# test_idx = np.arange(600)
# test_idx = T0[160:165]

print('test size %.3f sec\n' % ((test_idx[-1]-test_idx[0])*1e-3*dry_data.binwidth))

# Train on DRR sequences (inputs) and DRY targets!
_, y_train = dry_data[train_idx]    
X_train, _ = drr_data[train_idx]

# Test 
_, y_test = dry_data[test_idx]
X_test, y_test_drr = drr_data[test_idx]


trainset = TensorDataset( X_train, y_train )
testset = TensorDataset( X_test, y_test )

train_loader = DataLoader(trainset, batch_size=batch_size, shuffle=False)
test_loader = DataLoader(testset, batch_size=batch_size, shuffle=False)


# ! DEBUG
print('*** DEBUGING! ***')
X, y = next(iter(train_loader))
print('- X_train :', X_train.shape)
print('- y_train :', y_train.shape)
print('- X_test  :', X_test.shape)
print('- y_test  :', y_test.shape)



# %% #* Create RNN model
class RNN(nn.Module):
    def __init__(self, n_input, n_hidden, n_layers, 
        nonlinearity='relu', batch_first=True, dropout=0) -> None:
        super(RNN, self).__init__()
        self.n_input      = n_input         # dimensions of the input tensor
        self.n_hidden     = n_hidden        # num of hidden parameters
        self.n_layers     = n_layers        # num of recurrent RNN layers
        self.nonlinearity = nonlinearity    # {'tanh', 'relu'}
        self.batch_first  = batch_first     # input and output tensors are provided as (batch, seq, feature)
        self.dropout      = dropout 
        self.n_output     = n_input 

        # (num_layers * num_directions, batch_size, hidden_size)
        #self.hn = torch.randn(self.n_layers, batch_size, self.n_hidden).to(device)

        # RNN   # * hyper-parameter 
        self.rnn = nn.RNN(input_size=n_input, hidden_size=n_hidden, num_layers=n_layers,
           nonlinearity=nonlinearity, batch_first=batch_first, dropout=dropout)

        # GRU   # * hyper-parameter 
        #self.rnn = nn.GRU(input_size=n_input, hidden_size=n_hidden, num_layers=n_layers,
        #   batch_first=batch_first, dropout=dropout)

        self.fc = nn.Linear(n_hidden, self.n_output)

    def forward(self, x):
        # (num_layers * num_directions, batch_size, hidden_size)
        #self.h0 = torch.randn(self.n_layers, x.shape[0], self.n_hidden).to(device)
        self.h0 = torch.zeros(self.n_layers, x.shape[0], self.n_hidden).to(device)
    
        # ASSUMING batch_first == True...
        #   x input : (batch, seq_len, input_size)
        #   x output: (seq_len, batch, num_directions * hidden_size)
        x, hn = self.rnn(x, self.h0)
        x = self.fc( x[:,-1,:] )          
        return x


n_input = dry_data.n_units                      #* number of UNITs 
n_hidden = 100                                  #* hyper-parameter 
n_layers = 1                                    #* hyper-parameter 
nonlinearity = 'relu' # {'relu', 'tanh'}        #* hyper-parameter 
batch_first = True
dropout = 0
model = RNN(n_input = n_input, 
            n_hidden = n_hidden, 
            n_layers = n_layers,             
            nonlinearity = nonlinearity, 
            batch_first = True, 
            dropout = 0).to(device)


print(model)

# ! DEBUG
# print('*** DEBUGING! ***')
# X, y = next(iter(train_loader))
# yhat = model(X)
# print('- y.shape   : ', y.shape)
# print('- yhat.shape: ', yhat.shape)
# # X, y = dry_data[0]
# # yhat = model(X[None,:])     # add the BATCH dimension, as in (batch, seq_len, input_size)



# %%  # *** Training ***
epochs = 400

# Option 1
lr = 1e-3       # * hyper-parameter; learning rate
loss_function = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=lr)

# Option 2
# learning_rate = 0.05
# loss_function = nn.L1Loss()
# optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)


model.train()
loss = []

for ep in range(epochs):
    for k, (X, y) in enumerate(train_loader):
        optimizer.zero_grad()

        # ? Init the hidden states!!
         
        y_pred = model( X.to(device) )                
        loss_k = loss_function(y_pred.view(-1,1), y.view(-1,1).to(device))
        loss_k.backward()
        optimizer.step()


    loss.append(loss_k.item())
    if ep%25 == 1:
        print(f'epoch: {ep:4}, loss: {loss[-1]}')

print('- Finished training!')



# %% 
# *** Testing ***
model.eval()

y_est = []
y_dry = []

for k, (X, yk) in enumerate(test_loader):
    with torch.no_grad():   # ? DO I NEED no_grad() if I have eval() ? 
        yk_pred = model( X.to(device) )
        y_est.append(yk_pred.cpu().numpy())
        y_dry.append(yk.numpy())        

# 
y_est = np.concatenate(y_est, axis=0)
y_dry = np.concatenate(y_dry, axis=0)
# y_test_drr = y_test_drr.flatten()

assert np.all(np.equal(y_dry, y_test).numpy()), 'y_dry & y_test should be the same!'

CCdry = np.corrcoef(y_est.flatten(), y_dry.flatten())[0,1]
CCdrr = np.corrcoef(y_est.flatten(), y_test_drr.flatten())[0,1]

print(f'CCdry: {CCdry:8.3f}')
print(f'CCdrr: {CCdrr:8.3f}')

plt.figure(figsize=(14,8))
# plt.plot(np.c_[y_est, y_dry, y_drr])
plt.plot(y_dry[:500,0], label='$y_{dry}$')
plt.plot(y_est[:500,0], label='$\hat{y}$')
plt.plot(y_test_drr[:500,0], label='$y_{drr}$')
plt.legend()
plt.xlabel('Samples')
plt.ylabel('Amp.')
plt.title(f'$CC_{{dry-est}}: {CCdry:.3f}$, $CC_{{drr-est}}: {CCdrr:.3f}$')
plt.show()



# %%
# np.corrcoef(whitening(Hdrr[-30:]), whitening(y_drr), '.')[0,1]
# np.corrcoef(whitening(Hdry[-30:]), whitening(y_dry), '.')[0,1]

# %%
