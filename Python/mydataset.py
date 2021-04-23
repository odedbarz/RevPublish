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


import matplotlib.pyplot as plt
import os
from scipy.io import loadmat
import numpy as np

import torch
from torch.utils.data import Dataset
import pandas as pd


# %%
class dataset(Dataset):
    def __init__(self, seq_len, drr_idx, unit_numbers, fn, path='.', normalize=True, 
    cumsum=False) -> None:
        self.fn = fn        
        self.path = path

        assert os.path.isdir(self.path), 'No such directory!'
        full_path = os.path.join(self.path, self.fn)        
        assert os.path.isfile(full_path), 'No such file!'

        # Loads data from the MAT file
        dummy = loadmat(full_path, squeeze_me=True)

        self.drr_idx = drr_idx                                  # drr case to use
        self.drr_label = dummy['drr_labels'][self.drr_idx]      # drr label        
        self.unit_numbers = unit_numbers                          # unit number
        self.binwidth = dummy['binwidth']                       # (ms) binwidth

        # Get the DRR case and transform to torch tensor array
        M = torch.squeeze( torch.Tensor( dummy['data'] )[:,self.drr_idx,:] )

        # Extract desired units
        self.n_time, max_units = M.shape     # (time x number of units)    
        if np.isscalar(self.unit_numbers) and 0 <= self.unit_numbers:
            # Create a COLUMN vector of MUA/SU response
            self.response = M[:,self.unit_numbers][:,None]   
        elif 0==len(self.unit_numbers):            
            # Get all responses
            self.response = M   
            self.unit_numbers = (0, M.shape[-1])
        elif 2 == len(self.unit_numbers):
            self.response = M[:,self.unit_numbers[0]:self.unit_numbers[1]]
        else:
            raise Exception('unit_number must be of size 1 or 2 (scalar or matrix)!')
        
        self.cumsum = cumsum
        if self.cumsum:
            self.response = np.cumsum(self.response, axis=0)

        self.n_units = self.response.shape[-1]          
 
        # * NORMALIZE
        self.normalize = normalize
        if normalize:
            mu = self.response.mean(dim=0)
            std = self.response.std(dim=0)
            self.response = (self.response - mu)/std

        # * input of shape for RNNs, (seq_len, batch, input_size)
        self.seq_len = seq_len
        self.data_len = int(self.n_time-seq_len)
        self.X = torch.zeros(self.data_len, seq_len, self.n_units)      # a window of SEQ_LEN samples
        self.y = torch.zeros(self.data_len, self.n_units)               # next sample ("labels")
        for k in range(self.data_len):
            self.X[k,:] = self.response[k:k+seq_len,:]
            self.y[k] = self.response[k+seq_len,:] 

    def __repr__(self):
        '''
        A nice print of the object's fields
        '''
        str = 'Struct fields:\n'
        for key, value in self.__dict__.items():
            if isinstance(value, np.ndarray) or isinstance(value, torch.Tensor):
                new_var = value.shape
                value = ['shape: ', new_var]

            str += '- {0}: {1}\n'.format(key, value)
        return str
        
        # D = pd.DataFrame();
        # for key, value in self.__dict__.items():
        #     if isinstance(value, np.ndarray) or isinstance(value, torch.Tensor):
        #         value = value.shape
        #     print(key, value)
        #     D.append( pd.DataFrame([value], index=[key]) )
        
        # return D.to_string()

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx,:], self.y[idx]


# # %%
if __name__ == '__main__':
    print('Testing madataset.py...')

    # Example of loading only one unit response
    data_oneunit = dataset( seq_len = 20, drr_idx = 0, unit_numbers = 1, 
        fn = 'data_MUA-ONLY_(08-Jan-2021)_bw(5)_fbands(30)_spec(gammatone).mat',
        path = './')
    
    print(data_oneunit)


    # Example of loading all available unit responses
    data_allunits = dataset( seq_len = 20, drr_idx = 0, unit_numbers = [], 
        fn = 'data_MUA-ONLY_(08-Jan-2021)_bw(5)_fbands(30)_spec(gammatone).mat',
        path = './')
    
    print(data_allunits)
 

# %%

