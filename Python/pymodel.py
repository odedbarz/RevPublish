# %%
%reload_ext autoreload
%autoreload 2

from aux import whitening
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

unit_numbers = 0
dry_data = dataset( 
    seq_len = 20,    # binwidth * seq_len => duration in msec    # * hyper-parameter 
    drr_idx = 0, 
    unit_numbers = unit_numbers, 
    fn = fn, 
    path = path)

drr_data = dataset(
    seq_len = 20,    # binwidth * seq_len => duration in msec    # * hyper-parameter 
    drr_idx = 4, 
    unit_numbers = unit_numbers, 
    fn = fn, 
    path = path)

print(drr_data)


 
# %% Split the data into train/test sets
test_size = 0.2
batch_size = 64
shuffle = False

train_idx, test_idx = train_test_split(range(len(dry_data)), test_size=test_size, 
    shuffle=shuffle )  #, random_state=42)

# X_train, X_test, y_train, y_test = train_test_split(drr_data.X, dry_data.y, 
#     test_size=test_size, shuffle=shuffle ) #, random_state=42)

# # Compare with ...
# _, _, _, y_test_drr = train_test_split(np.empty_like(dry_data.X), drr_data.y, 
#     test_size=test_size, shuffle=shuffle) #, random_state=42)

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



# %% # * Create RNN model
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
        self.n_output     = n_output 

        # RNN   # * hyper-parameter 
        self.rnn = nn.RNN(input_size=n_input, hidden_size=n_hidden, num_layers=n_layers,
            nonlinearity=nonlinearity, batch_first=batch_first, dropout=dropout)

        # GRU   # * hyper-parameter 
        # self.rnn = nn.GRU(input_size=n_input, hidden_size=n_hidden, num_layers=n_layers,
        #    batch_first=batch_first, dropout=dropout)

        self.fc = nn.Linear(n_hidden, n_output)

    def forward(self, x):
       # (num_layers * num_directions, batch_size, hidden_size)
        self.h0 = torch.randn(self.n_layers, x.shape[0], self.n_hidden)

        # (batch, seq_len, input_size); for batch_first == True
        #x = x.view(x.shape[0], x.shape[1], 1)  

        # ASSUMING batch_first == True...
        #   x input : (batch, seq_len, input_size)
        #   x output: (seq_len, batch, num_directions * hidden_size)
        x, hn = self.rnn(x, self.h0)
        x = self.fc( x[:,-1,:] )          
        return x #[-1]


n_input = np.min((15, dry_data.n_units))
n_hidden = 100                                  # * hyper-parameter 
n_layers = 1                                    # * hyper-parameter 
n_output = 1
nonlinearity = 'relu' # {'relu', 'tanh'}
batch_first = True
dropout = 0
model = RNN(n_input = n_input, 
            n_hidden = n_hidden, 
            n_layers = n_layers, 
            n_output = n_output,
            nonlinearity = nonlinearity, 
            batch_first = True, 
            dropout = 0).to(device)


print(model)

# ! DEBUG
print('*** DEBUGING! ***')
X, y = next(iter(train_loader))
yhat = model(X)
print('- y.shape   : ', y.shape)
print('- yhat.shape: ', yhat.shape)
# X, y = dry_data[0]
# yhat = model(X[None,:])     # add the BATCH dimension, as in (batch, seq_len, input_size)



# %%  # *** Training ***
epochs = 400
lr = 1e-3       # learning rate

loss_function = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=lr)

model.train()
loss = []

for ep in range(epochs):
    for k, (X, y) in enumerate(train_loader):
        optimizer.zero_grad()

        # ? Init the hidden states!!
         
        y_pred = model( X.to(device) )
        # y_pred.squeeze_()
        loss_k = loss_function(y_pred, y)
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
        y_est.append(yk_pred.numpy())
        y_dry.append(yk.numpy())        

y_est = np.concatenate(y_est, axis=0).flatten()
y_dry = np.concatenate(y_dry, axis=0).flatten()
y_test_drr = y_test_drr.flatten()

assert np.all(np.equal(y_dry, y_test.flatten()).numpy()), 'y_dry & y_test should be the same!'

CCdry = np.corrcoef(y_est, y_dry)[0,1]
CCdrr = np.corrcoef(y_est, y_test_drr)[0,1]

print(f'CCdry: {CCdry:8.3f}')
print(f'CCdrr: {CCdrr:8.3f}')

# plt.plot(np.c_[y_est, y_dry, y_drr])
plt.plot(y_dry, label='$y_{dry}$')
plt.plot(y_est, label='$\hat{y}$')
plt.plot(y_test_drr, label='$y_{drr}$')
plt.legend()
plt.xlabel('Samples')
plt.ylabel('Amp.')
plt.title(f'$CC_{{dry-est}}: {CCdry:.3f}$, $CC_{{drr-est}}: {CCdrr:.3f}$')




# %%
# np.corrcoef(whitening(Hdrr[-30:]), whitening(y_drr), '.')[0,1]
# np.corrcoef(whitening(Hdry[-30:]), whitening(y_dry), '.')[0,1]

