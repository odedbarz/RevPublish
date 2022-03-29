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
n_wav_files = 4620;     % (# of files in train: 4621) number of wav files to use
sr          = 16000;    % wav sample rate

spectrogram_type = 'gammatone';      % {['matlab'], 'stft', 'multitaper', 'gammatone'}
f_scale     = 'erb';        % {['lin'], 'log', 'erb'}
n_freq      = 64, '**!!**'  % (1x1) # of bins along the frequency domain of the spectrogram
binwidth    = 5             % (ms) binwidth of the resulted spectrogram 
win_size_ms = nan           % (ms) temporal window size over which to calc the spectrogram; 
                            %      'gammatone' filterbanks do not use it!
lowfreq     = 250;          % (Hz)
highfreq    = 8000;         % (Hz) 8900 Hz so that to get as high as possible to near 8k Hz with the 
nw          = [];           % applies only for SPECTROGRAM_TYPE = 'multitaper'
db_floor    = -80;

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


%%
too_short_files = 0;

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
        y_left = y_conv{n}(:,1);    % LEFT RIR        
        [Sft, spec_st] = spec.spectrogram(y_left, fs, ...
            'n_bands', n_freq,...
            'lowfreq', lowfreq,...
            'highfreq', highfreq,...
            'overlap_ratio', nan,...
            'binwidth', binwidth,...
            'win_size_ms', win_size_ms, ...
            'nw', nw,...                only for spectrogram_type== MULTITAPER
            'f_scale', f_scale,...
            'db_floor', db_floor, ...  % (dB)
            'duration_ms', duration_ms,...
            'method', spectrogram_type, ...
            'fignum', []...
         ); 
     
        
        y_right = y_conv{n}(:,2);    % RIGHT RIR    
        [Sft_right, spec_st] = spec.spectrogram(y_right, fs, ...
            'n_bands', n_freq,...
            'lowfreq', lowfreq,...
            'highfreq', highfreq,...
            'overlap_ratio', nan,...
            'binwidth', binwidth,...
            'win_size_ms', win_size_ms, ...
            'nw', nw,...                only for spectrogram_type== MULTITAPER
            'f_scale', f_scale,...
            'db_floor', db_floor, ...  % (dB)
            'duration_ms', duration_ms,...
            'method', spectrogram_type, ...
            'fignum', []...
         );       
     
    
%         if size(Sft, 2) < MAX_SAMPLES
%             too_short_files = too_short_files + 1;
%             fprintf('- too_short_files: %g (k=%d)\n', too_short_files, k);
%             break;
        %else
        %   Sft = Sft(:,1:MAX_SAMPLES);
%         end        
     
        save_to_label = fullfile( save_to_path, drr_labels{n} );
        if ~isfolder( save_to_label )
            mkdir( save_to_label );
        end     
        fn = sprintf('%04d_Sft_drr(%s)_freq(%d).mat', k,...
            drr_labels{n}, n_freq);
        fn = fullfile( save_to_label, fn );
        save(fn, 'Sft', 'Sft_right');
        
    end
 
 
end



%% Save the parameters
fn = sprintf('datainfo_(%s)_freq(%d).mat', date, n_freq);
save( fullfile(save_to_path, fn), ...
    'spec_st', 'Trir');







