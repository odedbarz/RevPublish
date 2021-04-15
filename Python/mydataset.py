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
    def __init__(self, drr_idx, unit_number, fn, path='.') -> None:
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
        self.v = M[:,self.unit_number]
        self.binwidth = dummy['binwidth']       # (ms) binwidth

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
        return self.v[idx]


# # %%
if __name__ == '__main__':
    print('Testing madataset.py...')

    data = dataset( drr_idx = 0, unit_number = 1, 
        fn = 'data_MUA-ONLY_(08-Jan-2021)_bw(5)_fbands(30)_spec(gammatone).mat',
        path = './')
    
    print(data)

    print('data[:5]')
    print( torch.stack((data[:5], data.v[:5]), dim=1) )
 
# %%
