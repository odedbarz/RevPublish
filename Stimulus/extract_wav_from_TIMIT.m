% 
% extract_wav_from_TIMIT.m
% 

clc

timit_path  = 'D:\DataBases\TIMIT\';
path_subs   = 'TIMIT';
save_2_path = 'D:\DataBases\TIMIT\wav_files\';
n_files     = 2500;    % # of files to load from each gender

fprintf('--> Starting to extract WAV files from the TIMIT corpus...\n');
lookfor_wavs_in_subfolders(timit_path, path_subs, save_2_path, n_files);
fprintf('--> Finished!\n');





