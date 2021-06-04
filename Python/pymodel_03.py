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
from matplotlib.pyplot import *
import matplotlib
matplotlib.style.use('seaborn')

import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

from sklearn.model_selection import train_test_split

from mydataset import dataset

# set device
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f'device: {device}')



# %% Load from MAT file
fn = 'data_MUA-ONLY_(08-Jan-2021)_bw(5)_fbands(30)_spec(gammatone).mat'     # MUA
# fn = 'data_SU-ONLY_(08-Jan-2021)_bw(5)_fbands(30)_spec(gammatone).mat'      # SU 
# fn = 'data_SU-Hann(25)_(23-Apr-2021)_bw(5)_fbands(30)_spec(gammatone).mat'
path = './'

seq_len = 50                 # samples   #* hyper-parameter 
unit_numbers = 18                    #* NUMBER of UNITS for analysis 
normalize = True

dry_data = dataset( 
    seq_len = seq_len,      # binwidth * seq_len => duration in msec    
    drr_idx = 0, 
    unit_numbers = unit_numbers, 
    fn = fn, 
    path = path,
    normalize=normalize, 
    cumsum=False)

# fn_strf = 'STRF_MUA-ONLY_bw(5)_fbands(30)_spec(gammatone).mat'
drr_data = dataset(
    seq_len = seq_len,      # binwidth * seq_len => duration in msec    
    drr_idx = 4, 
    unit_numbers = unit_numbers, 
    fn = fn, 
    path = path, 
    normalize=normalize,
    cumsum=False)

print(drr_data)


 
# %% Split the data into train/test sets
test_size = 0.1
batch_size = 64*2
shuffle = True

train_idx, test_idx = train_test_split(range(len(dry_data)), test_size=test_size, 
    shuffle=shuffle )  #, random_state=42)

print('test size %.3f sec\n' % ((test_idx[-1]-test_idx[0])*1e-3*dry_data.binwidth))

# Train on DRR sequences (inputs) and DRY targets!
_, y_train = dry_data[train_idx]    
X_train, _ = drr_data[train_idx]
X_train = X_train.permute(0,2,1)      # (batch, units, seq) 

# Test 
_, y_test = dry_data[test_idx]
X_test, y_test_drr = drr_data[test_idx]
X_test = X_test.permute(0,2,1)      # (batch, units, seq)

trainset = TensorDataset(X_train, y_train)
testset = TensorDataset(X_test, y_test)

train_loader = DataLoader(trainset, batch_size=batch_size, shuffle=False)
test_loader = DataLoader(testset, batch_size=batch_size, shuffle=False)


# ! DEBUG
print('*** DEBUGING! ***')
X, y = next(iter(train_loader))
print('- X_train :', X_train.shape)
print('- y_train :', y_train.shape)
print('- X_test  :', X_test.shape)
print('- y_test  :', y_test.shape)



# %%
class Conv1d(nn.Module):
    def __init__(self, n_input=1, kernel_size=3) -> None:
        super().__init__()
        self.n_input = n_input
        self.n_output = n_input
        self.kernel_size = kernel_size
        self.p = 0.1

        self.conv1 = nn.Conv1d(self.n_input, 5, self.kernel_size, padding=0, bias=True)
        self.relu = nn.ReLU()
        self.maxpool1d = nn.MaxPool1d(kernel_size=3, stride=2, padding=1)
        self.batchnorm = nn.BatchNorm1d(num_features=5)
        self.layer1 = nn.Sequential(self.conv1, self.relu, self.maxpool1d, self.batchnorm)
        
        self.conv2 = nn.Conv1d(5, 10, 3, padding=2, bias=True)
        # self.relu = nn.ReLU()
        # self.maxpool1d = nn.MaxPool1d(kernel_size=3, stride=2, padding=1)
        # self.batchnorm = nn.BatchNorm1d(num_features=5)
        self.layer2 = nn.Sequential(self.conv2, self.relu, self.maxpool1d, 
            nn.BatchNorm1d(num_features=10))
 
        self.dropout = nn.Dropout(p=self.p)
        
        self.fc_layer1 = nn.Sequential(
            nn.Linear(120, 1), self.dropout
        )


    def forward(self, x):
        # input size: (batch_size, channels, seq_length)
        
        # #! DEBUG
        # x = self.layer1(x)
        # x = self.layer2(x)
        # print('- x.shape: ', x.shape)
        # x = self.conv2(x)
        # print('- x.shape: ', x.shape)
        # x = self.relu(x)
        # x = self.maxpool1d(x)
        # print('- x.shape: ', x.shape)
        # x = self.batchnorm(x)
        # print('- x.shape: ', x.shape)

        #self.conv.padding = 0
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.fc_layer1(x.view(-1,x.shape[1]*x.shape[2]))
        return x

# model = Conv1d(n_input=dry_data.n_units, kernel_size=7)
# print(model)


# ! DEBUG
X, y = next(iter(train_loader))
y_ = model(X)
print('- X : ', X.shape)
print('- y : ', y.shape)
print('- y_: ', y_.shape)


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

        self.hn = []

        # (num_layers * num_directions, batch_size, hidden_size)
        #self.hn = torch.randn(self.n_layers, batch_size, self.n_hidden).to(device)

        # RNN   # * hyper-parameter 
        #self.rnn = nn.RNN(input_size=n_input, hidden_size=n_hidden, num_layers=n_layers,
        #   nonlinearity=nonlinearity, batch_first=batch_first, dropout=dropout)

        # GRU   # * hyper-parameter 
        self.rnn = nn.GRU(input_size=n_input, hidden_size=n_hidden, num_layers=n_layers,
          batch_first=batch_first, dropout=dropout)

        #self.fc = nn.Linear(n_hidden, self.n_output)
        self.fc1 = nn.Linear(n_hidden, self.n_output)

    def init_hidden(self, x):
         #self.h0 = torch.randn(self.n_layers, x.shape[0], self.n_hidden).to(device)
        self.hn = torch.zeros(self.n_layers, x.shape[0], self.n_hidden).to(device)       
        return None

    def forward(self, x):
        x = x.permute(0,2,1)

        # (num_layers * num_directions, batch_size, hidden_size)
        #self.h0 = torch.randn(self.n_layers, x.shape[0], self.n_hidden).to(device)
        #self.h0 = torch.zeros(self.n_layers, x.shape[0], self.n_hidden).to(device)

        # ASSUMING batch_first == True...
        #   x input : (batch, seq_len, input_size)
        #   x output: (seq_len, batch, num_directions * hidden_size)
        x, self.hn = self.rnn(x, self.hn)
        #x = self.fc( x[:,-1,:] )   
        x = self.fc1(x)   
        #x = self.fc2(x.squeeze()) 
        x = x[:,-1]
        return x

n_input = dry_data.n_units                      #* number of UNITs 
n_hidden = 100                                  #* hyper-parameter 
n_layers = 1                                    #* hyper-parameter 
nonlinearity = 'relu' # {'relu', 'tanh'}        #* hyper-parameter 
batch_first = True
dropout = 0.3
model = RNN(n_input = n_input, 
            n_hidden = n_hidden, 
            n_layers = n_layers,             
            nonlinearity = nonlinearity, 
            batch_first = batch_first, 
            dropout = 0).to(device)

print(model)

# ! DEBUG
print('*** DEBUGING! ***')
X, y = next(iter(train_loader))
model.init_hidden(X)
yhat = model(X)
print('- y.shape   : ', y.shape)
print('- yhat.shape: ', yhat.shape)
# X, y = dry_data[0]
# yhat = model(X[None,:])     # add the BATCH dimension, as in (batch, seq_len, input_size)

#! DEBUG
# data = loadmat('../_data/Analysis/analyzed_librosa.mat', squeeze_me=True)



# %%  # *** Training ***
epochs = 60

# - OPTIMIZERs
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=0)
# optimizer = torch.optim.Adadelta(model.parameters())
# optimizer = torch.optim.SGD(model.parameters(), lr=0.001)


# - LOSS FUNCTIONs
loss_function = nn.MSELoss()
# loss_function = nn.L1Loss()
# def loss_function(a, b):
#     av = a.flatten() - a.mean()
#     bv = b.flatten() - b.mean()
#     cov = torch.dot(av, bv)/(torch.norm(av)*torch.norm(bv))
#     return 1.0 - cov


model.train()
loss = []

for ep in range(epochs):
    #model.init_hidden()
    for k, (X, y) in enumerate(train_loader):
        model.init_hidden(X)
        optimizer.zero_grad()
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
    with torch.no_grad():   
        model.init_hidden(X)
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
CCdry2drr = np.corrcoef(y_dry.flatten(), y_test_drr.flatten())[0,1]

print(f'CCdry    : {CCdry:8.3f}')
print(f'CCdrr    : {CCdrr:8.3f}')
print(f'CCdry2drr: {CCdry2drr:8.3f}')

plt.figure(figsize=(14,8))
# plt.plot(np.c_[y_est, y_dry, y_drr])
plt.plot(y_dry[:500,0], label='$y_{dry}$')
plt.plot(y_est[:500,0], label='$\hat{y}$')
plt.plot(y_test_drr[:500,0], label='$y_{drr}$')
plt.legend()
plt.xlabel('Samples')
plt.ylabel('Amp.')
plt.title(f'$CC_{{dry-est}}: {CCdry:.3f}$, $CC_{{drr-est}}: {CCdrr:.3f}$, $CC_{{dry-drr}}: {CCdry2drr:.3f}$')
plt.show()

# plt.figure(figsize=(14,8))
# plot(whitening( np.c_[y_dry[:500,0], y_est[:500,0], y_test_drr[:500,0]])) 
# plt.show()

# %%
