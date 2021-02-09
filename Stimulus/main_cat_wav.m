%
% main_cat_wav.m
%
% Description:
% * Adds (aggregates) together all wav files into one big file;
% * Each WAV file is seperated by silence.
%

% Fs = 100e3; 	% (Hz)
Fs = 100e3; 	% (Hz)

% The new file name to save
fn_dry = 'TIMIT_dry';

% Define paths for the WAV files
% path_females = './TIMIT_Females/';
% path_males   = './TIMIT_Males/';
path_females = './TIMIT_Females_Fs(100kHz)/';
path_males   = './TIMIT_Males_Fs(100kHz)/';

% Concatenate the signals
[yv, seg_time_samples] = cat_wav(...
    'Fs_base', Fs,...
    'path_females', path_females, ...
    'path_males', path_males, ...
    'fn_dry', fn_dry ...    % SAVE the concatenated file
);

save([fn_dry, '_meta.mat'], 'Fs', 'seg_time_samples');














