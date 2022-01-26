%
% main_aggregate_SU_data.m
%
% Description:
% This script run over all measurements and aggregate all the PSTH data into
% one big matrix H. 
%
%   H: (time-samples, DRR-cases, # of measurements)
%
% The data is then saved into a MAT file.
%

clc
fignum = 11;
verbose = 1;

setup_environment('../');


drr     = get_DRR_list_and_indices;
n_drr 	= drr.n_drr;  


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('et_data_files');



%% Load the measurement's table
if ~exist('tbl_main', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Row is a duplicate variable name used by MATLAB
    tbl_impale.Properties.VariableNames{'Row'} = 'neuron';    
end
    


%% Load stimuli & measurement data
fprintf('--> Load stimuli & spectrograms\n');
spectrogram_type = 'gammatone';      % {['matlab'], 'stft', 'multitaper', 'gammatone', 'meddis', 'carney'}
f_scale     = 'erb';	% {['lin'], 'log', 'erb'}
n_bands     = 3*30      % (1x1) # of bins along the frequency domain of the spectrogram
binwidth    = 1         % (ms) binwidth of the resulted spectrogram 
amp_to_db   = true
win_size_ms = nan       % (ms) temporal window size over which to calc the spectrogram; 
                        %      'gammatone' filterbanks do not use it!
lowfreq     = 125;      % (Hz)
highfreq    = 8000;   % (Hz) 8900 Hz so that to get as high as possible to near 8k Hz with the 
nw          = [];       % applies only for SPECTROGRAM_TYPE = 'multitaper'
duration_to_load = 36;  % neuron #1 (#115) is for stimulus of 36 sec (40 sec)

spectral_diff     = 0;       % (logical) perform derivative (DIFF) along the frequency domain
apply_sync_filter = false
calc_envelope     = true

% STIM_LIST & SPEC_LIST contains all stimuli duration (i.e., the 36 s &
% 40 s stimuli)
[stim_st, spec_st] = load.stimulus_and_spectrogram(tbl_impale, ...
    'spectrogram_type', spectrogram_type, ...
    'duration_to_load', duration_to_load, ...
    'binwidth', binwidth,...
    'lowfreq', lowfreq, ...         % (Hz)
    'highfreq', highfreq, ...       % (Hz)
    'win_size_ms', win_size_ms, ...
    'n_bands', n_bands, ...
    'f_scale', f_scale, ... {'lin', 'log'}
    'nw', nw, ...          (default: 1.4) only for spectrogram_type == MULTITAPER
    'spectral_diff', spectral_diff, ...
    'apply_sync_filter', apply_sync_filter, ...
    'calc_envelope', calc_envelope, ...
    'fignum', [] ...
    );            

fs = 1/(1e-3*binwidth);         % (Hz)

% ALways load stimuli with the same duration
'Loads stimuli with the SAME duration'
duration_sec= 36;   % (sec) stimulus duration to use 
duration_ms = units.sec2ms( duration_sec );

% Plot spectrogram
spec.plot_spectrogram(spec_st.t, spec_st.f, spec_st.Sft{3});
title(sprintf('Spectrogram (DRY)'));

neuron_idx  = tbl_impale.SPK & (tbl_impale.duration_sec == duration_sec);
neuron_list = find(neuron_idx);
tbl_SU      = tbl_impale(neuron_idx, :);
fprintf('--> NOTE: tbl_SU contains stimuli of %g sec ONLY!\n', unique(tbl_SU.duration_sec));


 
%% Load PSTH responses
S_list = load.response( tbl_SU, duration_sec );
 
if verbose
    aux.cprintf('UnterminatedStrings', '\n    Stimulus Parameters:\n');    
    aux.cprintf('UnterminatedStrings', '--> Sampling rate (fs): %g kHz\n', 1e-3*stim_st.fs);
    aux.cprintf('UnterminatedStrings', '--> Selected DURATION : %g (sec)\n', 1e-3*duration_sec);    
    aux.cprintf('UnterminatedStrings', '\n    Spectrogram Parameters:\n');
    aux.cprintf('UnterminatedStrings', '--> spectrogram_type: %s\n', spectrogram_type);    
    aux.cprintf('UnterminatedStrings', '--> binwidth        : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands         : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms     : %g ms\n', win_size_ms);    
    aux.cprintf('UnterminatedStrings', '--> f_scale         : %s\n', f_scale);
    aux.cprintf('UnterminatedStrings', '--> lowfreq         : %f Hz\n', lowfreq);
    aux.cprintf('UnterminatedStrings', '--> highfreq        : %f Hz\n', highfreq);
end




%% Load RESPONSE measureemnts
aux.vprint(verbose, '\n--> Load response H data...\n');

pst_type = {'psth'};
% pst_type = {'pstw', win_size_ms};    % {'pstw', WIN_SIZE_MS}
H = calc_PSTHs(S_list, binwidth, 'pst_type', pst_type);

aux.vprint(verbose, '--> Finished\n');



%% Get all rows (measurements) that have full session (all DRR conditions)
H_labels = zeros(size(H,3), size(H,2));
for unit = 1:size(H,3)     % # of measurement
    for ndrr = 1:size(H,2)     % DRR case
        H_labels(unit, ndrr) = any(~isnan( H(:,ndrr,unit) ));
    end
end

valid_neuron_idx = n_drr <= sum(H_labels(:,1:n_drr),2);
tbl_SU = tbl_SU(valid_neuron_idx, :);
H      = H(:,:,valid_neuron_idx);
S_list = S_list(valid_neuron_idx);




%% Save the data
'SAVE the analysis!'
fn.path = load.path_to_data('data');
fn.file = sprintf('data_SU_(%s)_bw(%g)_fbands(%d)_win(%g)ms_spec(%s).mat', ...
    date, binwidth, n_bands, win_size_ms, spectrogram_type);
fn.save = fullfile(fn.path, fn.file);

fprintf('\nSaving data at:\n');
disp(fn)
save(fn.save, 'H', 'tbl_SU', 'spec_st', 'stim_st', 'S_list');







