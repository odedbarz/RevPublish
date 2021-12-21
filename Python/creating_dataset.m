%
% main_aggregate_MUA_data.m
%
% Description:
% This script run over all measurements and aggregate all the MUA data into
% one big matrix H. 
%
%   H: (time-samples, DRR-cases, # of measurements)
%
% The data is then saved into a MAT file.
%

clc
fignum = 11;
verbose = 1;

addpath('..');
addpath(fullfile('..', 'Stimulus'));    % for ConvRIR
setup_environment('..');



%% Get a list of ALL available wav files
timit_root_path = fullfile('..', '..', 'PyTorch', '_datasets', 'TIMIT', 'train');
wav_files = dir( fullfile(timit_root_path, '**', '*.wav') );



%%
n_wav_files = 2000;     % number of wav files to use
sr          = 16000;    % wav sample rate

spectrogram_type = 'gammatone';      % {['matlab'], 'stft', 'multitaper', 'gammatone'}
f_scale     = 'erb';        % {['lin'], 'log', 'erb'}
n_freq      = 224, '**!!**'  % (1x1) # of bins along the frequency domain of the spectrogram
binwidth    = 1             % (ms) binwidth of the resulted spectrogram 
win_size_ms = nan           % (ms) temporal window size over which to calc the spectrogram; 
                            %      'gammatone' filterbanks do not use it!
lowfreq     = 100;          % (Hz)
highfreq    = 8100;         % (Hz) 8900 Hz so that to get as high as possible to near 8k Hz with the 
nw          = [];           % applies only for SPECTROGRAM_TYPE = 'multitaper'


%% Loads the RIR responses
dummy = load( fullfile( load.path_to_data('Stimulus'), 'Spch_(36)sec/spch_36_metadata_new.mat') );
Trir = dummy.Trir
sort_by_drr = [1,2,5,3,6];      % sorted from DRY to most reverberant (-8.2 dB) RIR 
Trir = Trir(sort_by_drr, :);

drr = get_DRR_list_and_indices;
drr_labels = drr.labels([3, 2, 5, 1, 4]);

drr_to_use = [1, 5];
Trir = Trir(drr_to_use, :);
drr_labels = drr_labels(drr_to_use);

disp( Trir );



%%
save_to_path = sprintf('./_datasets/_data_freq(%d)', n_freq);
if ~isfolder( save_to_path )
    mkdir( save_to_path );
end


for k = 1:n_wav_files
    fn = fullfile( wav_files(k).folder, wav_files(k).name );
    assert(~isempty(dir(fn)), sprintf('ERROR: can''t find the file:\n\t==>> %s',fn));    
    
    [yk, fs] = audioread(fn);
    assert(fs==sr, sprintf('ERROR: There is a mismatch with the loaded sampling rate!'));
    file_info = audioinfo(fn);
    duration_ms = round( file_info.Duration * 1e3);     %(ms)
    
    % RIR
    y_conv = ConvRIR(yk, Trir.rir);
    y_conv = y_conv(:,1);   % take the *LEFT* size
    
    for n = 1:length(y_conv)
        yn = y_conv{n}(:,1);    % use the LEFT RIR only!

        [Sft, spec_st] = spec.spectrogram(yn, fs, ...
            'n_bands', n_freq,...
            'lowfreq', lowfreq,...
            'highfreq', highfreq,...
            'overlap_ratio', nan,...
            'binwidth', binwidth,...
            'win_size_ms', win_size_ms, ...
            'nw', nw,...                only for spectrogram_type== MULTITAPER
            'f_scale', f_scale,...
            'db_floor', -100, ...  % (dB)
            'duration_ms', duration_ms,...
            'method', spectrogram_type, ...
            'fignum', []...
         ); 

        save_to_label = fullfile( save_to_path, drr_labels{n} );
        if ~isfolder( save_to_label )
            mkdir( save_to_label );
        end     
        fn = sprintf('%03d_Sft_drr(%s)_freq(%d).mat', k, drr_labels{n}, n_freq);
        fn = fullfile( save_to_label, fn );
        %save(fn, 'Sft');
        
    end
 
 
end


fn = sprintf('datainfo_(%s)_freq(%d).mat', date, n_freq);
save( fullfile(save_to_path, fn), ...
    'spec_st', 'Trir');







