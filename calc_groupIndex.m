function groupIndex = calc_groupIndex(len_stim, binwidth, fs_stim, fn_data)
%
%   function calc_groupIndex(len_stim, binwidth, fs_stim, fn_data)
%
% Input: 
%   len_stim: (1x1) length of the concatenated input vector.
%   binwidth: (1x1, Hz) sampling rate of the spectrogram
%   fs_stim : (1x1, Hz) sampling rate of the stimuli
%   fn_data : (str) filename for the table that holds all info about the
%             TIMIT WAV files.
%
% Output:
%   groupIndex: (1 x len_stim) a vector of indices; each indice belongs to
%               another chunk of a TIMIT wav file.
%
% Description:
%   Assign an index to each chunk of the TIMIT stimuli that compose the
%   main concatenated stimulus used in the measurements.
%
%

if 3 > nargin
    fs_stim = 16e3;     % (Hz) teh default sampling rate of the TIMIT WAV files.
end

% Default filename
if 4 >nargin
    fn_data = [load.path_to_data('project1_new'), '_spch_metadata.mat'];
end


%%
% Load the TIMIT WAV's data table 
load(fn_data, 'tbl_data');

% Number of total samples
t = (0:len_stim)'*(1e-3*binwidth);

% # of TIMIT wav files used to create the (concatenated) whole stimulus
n_wavs = size(tbl_data,1);

% Split the (concatenated) stimulus into the TIMIT segments
groupIndex = nan(1, len_stim);
idx0 = 1;   % (smp)
for kk = 1:n_wavs
    % Translate from the stimulus' to the spectrogram's sampling rate
    samples_k = tbl_data.wav_end(kk);
    [~, idx1] = min(abs(samples_k/fs_stim - t));
    
    %idx1 = tbl_data.wav_end(kk);
    groupIndex(idx0:idx1) = kk;  
    
    idx0 = 1 + idx1;
end

% Add the last samples of 'silence'
groupIndex(idx0:end) = kk;






