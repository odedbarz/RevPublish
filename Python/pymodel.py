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

unit_number = 1
Hdry = dataset( drr_idx = 0, unit_number = unit_number, fn = fn, path = path)
Hdrr = dataset( drr_idx = 4, unit_number = unit_number, fn = fn, path = path)

# * NORMALIZE
mu = Hdry[:].mean()
std = Hdry[:].std()
Hdry.v = (Hdry.v - mu)/std
Hdrr.v = (Hdrr.v - mu)/std

print(Hdrr)

# Pack sequences and labels (next-sequence) into one list
seq_len = 20    # binwidth * seq_len => duration in msec    # * hyper-parameter 
data_length = int(Hdry.n_time-seq_len)
in_seq = torch.zeros(data_length, seq_len)      # in sequence
labels = torch.zeros(data_length)               # out sequence
labels_drr = torch.zeros(data_length)               # out sequence
for ii in range(data_length):
    in_seq[ii,:] = Hdrr[ii:ii+seq_len]
    labels[ii] = Hdry[ii+seq_len]
    labels_drr[ii] = Hdrr[ii+seq_len]
 
# Split the data into train/test sets
test_size = 0.2
X_train, X_test, y_train, y_test = train_test_split(in_seq, labels, test_size=test_size) #, random_state=42)
_, _, _, y_test0 = train_test_split(in_seq, labels_drr, test_size=test_size) #, random_state=42)

trainset = TensorDataset( X_train, y_train )
testset = TensorDataset( X_test, y_test )

batch_size = 64
shuffle = False
train_loader = DataLoader(trainset, batch_size=batch_size, shuffle=shuffle)
test_loader = DataLoader(testset, batch_size=batch_size, shuffle=shuffle)



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

        # (batch, seq_len, input_size); batch_first == True
        x = x.view(x.shape[0], x.shape[1], 1)       #! batch_size == len(x) == train_window

        x, hn = self.rnn(x, self.h0)
        x = self.fc( x[:,-1,:] )          
        return x #[-1]


# %% Create a model
n_input = 1
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
# X, y = next(iter(train_loader))
# yhat = model(X)  #! *TEST*
# plt.plot(np.c_[y.numpy(), yhat.detach().numpy().squeeze()])



# %% 
# *** Training ***
epochs = 400
lr = 1e-3       # learning rate

loss_function = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=lr)

model.train()
loss = []

for ep in range(epochs):
    for k, (seq, y) in enumerate(train_loader):
        optimizer.zero_grad()

        #! Init the hidden states!!
         
        y_pred = model( seq.to(device) )
        y_pred.squeeze_()
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
for k, (seq, yk) in enumerate(test_loader):
    with torch.no_grad():   # ? DO I NEED no_grad() if I have eval() ? 
        yk_pred = model( seq.to(device) )
        y_est.append(yk_pred.numpy())
        y_dry.append(yk.numpy())

idx = -1
y_est = whitening(y_est[idx])
y_dry = whitening(y_dry[idx])
# y_drr = whitening(labels_drr[-len(y_est):])
y_drr = whitening(y_test0[-len(y_est):])

CCdry = np.corrcoef(y_est, y_dry)[0,1]
CCdrr = np.corrcoef(y_est, y_drr)[0,1]

print(f'CCdry: {CCdry:8.3f}')
print(f'CCdrr: {CCdrr:8.3f}')

plt.plot(np.c_[y_est, y_dry, y_drr])
plt.plot(y_dry, label='$y_{dry}$')
plt.plot(y_est, label='$\hat{y}$')
plt.plot(y_drr, label='$y_{drr}$')
plt.legend()
plt.xlabel('Samples')
plt.ylabel('Amp.')
plt.title(f'$CC_{{dry-est}}: {CCdry:.3f}$, $CC_{{drr-est}}: {CCdrr:.3f}$')


# %%
# np.corrcoef(whitening(Hdrr[-30:]), whitening(y_drr), '.')[0,1]
# np.corrcoef(whitening(Hdry[-30:]), whitening(y_dry), '.')[0,1]



# %%
