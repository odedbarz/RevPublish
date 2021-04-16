# %%
# %reload_ext autoreload
# %autoreload 2

import matplotlib.pyplot as plt
import os
from scipy.io import loadmat
import numpy

import torch
from torch.utils.data import Dataset
import pandas as pd


# %%

class dataset(Dataset):
    def __init__(self, seq_len, drr_idx, unit_number, fn, path='.', normalize=True) -> None:
        self.fn = fn        
        self.path = path

        assert os.path.isdir(self.path), 'No such directory!'
        full_path = os.path.join(self.path, self.fn)        
        assert os.path.isfile(full_path), 'No such file!'

        # Loads data from the MAT file
        dummy = loadmat(full_path, squeeze_me=True)

        self.drr_idx = drr_idx                                  # drr case to use
        self.drr_label = dummy['drr_labels'][self.drr_idx]      # drr label        
        self.unit_number = unit_number                          # unit number

        # Get the DRR case and transform to torch tensor array
        M = torch.squeeze( torch.Tensor( dummy['mua'] )[:,self.drr_idx,:] )

        self.n_time, _ = M.shape                # length of measurement
        self.response = M[:,self.unit_number]          # vector of MUA/SU response
        self.binwidth = dummy['binwidth']       # (ms) binwidth

        # * NORMALIZE
        self.normalize = normalize
        if normalize:
            mu = self.response.mean()
            std = self.response.std()
            self.response = (self.response - mu)/std

        self.seq_len = seq_len
        self.data_len = int(self.n_time-seq_len)
        self.X = torch.zeros(self.data_len, seq_len)      # a window of SEQ_LEN samples
        self.y = torch.zeros(self.data_len)               # next sample ("labels")
        for k in range(self.data_len):
            self.X[k,:] = self.response[k:k+seq_len]
            self.y[k] = self.response[k+seq_len] 

    def __repr__(self):
        '''
        A nice print of the object's fields
        '''
        str = 'Struct fields:\n'
        for key, value in self.__dict__.items():
            if isinstance(value, numpy.ndarray) or isinstance(value, torch.Tensor):
                new_var = value.shape
                value = ['shape: ', new_var]

            str += '- {0}: {1}\n'.format(key, value)
        return str
        
        # D = pd.DataFrame();
        # for key, value in self.__dict__.items():
        #     if isinstance(value, numpy.ndarray) or isinstance(value, torch.Tensor):
        #         value = value.shape
        #     print(key, value)
        #     D.append( pd.DataFrame([value], index=[key]) )
        
        # return D.to_string()

    def __len__(self):
        return (self.n_time)

    def __getitem__(self, idx):
        return self.X[idx,:], self.y[idx]


# # %%
if __name__ == '__main__':
    print('Testing madataset.py...')

    data = dataset( seq_len = 20, drr_idx = 0, unit_number = 1, 
        fn = 'data_MUA-ONLY_(08-Jan-2021)_bw(5)_fbands(30)_spec(gammatone).mat',
        path = './')
    
    print(data)


