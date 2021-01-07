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



%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the measurement's table
if ~exist('tbl_main', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
end
    


%% Load stimuli & measurement data
fprintf('--> Load stimuli & spectrograms\n');
spectrogram_type = 'gammatone';      % {['matlab'], 'stft', 'multitaper', 'gammatone'}
f_scale     = 'log';	% {['lin'], 'log', 'erb'}
n_bands     = 30        % (1x1) # of bins along the frequency domain of the spectrogram
binwidth    = 5         % (ms) binwidth of the resulted spectrogram 
win_size_ms = nan       % (ms) temporal window size over which to calc the spectrogram; 
                        %      'gammatone' filterbanks do not use it!
lowfreq     = 250;      % (Hz)
highfreq    = 8840;     % (Hz) %8800;     
nw          = [];       % applies only for SPECTROGRAM_TYPE = 'multitaper'
neurons     = 1;        % neuron #1 (#115) is for stimulus of 36 sec (40 sec)
spectral_diff= 0;       % (logical) perform derivative (DIFF) along the frequency domain
hpf_pole    = nan;

% STIM_LIST & SPEC_LIST contains all stimuli duration (i.e., the 36 s &
% 40 s stimuli)
[stim_st, spec_st] = load.stimulus_and_spectrogram(tbl_impale, ...
    'spectrogram_type', spectrogram_type, ...
    'neurons', neurons, ...
    'binwidth', binwidth,...
    'lowfreq', lowfreq, ...         % (Hz)
    'highfreq', highfreq, ...       % (Hz)
    'win_size_ms', win_size_ms, ...
    'n_bands', n_bands, ...
    'f_scale', f_scale, ... {'lin', 'log'}
    'nw', nw, ...          (default: 1.4) only for spectrogram_type == MULTITAPER
    'spectral_diff', spectral_diff, ...
    'hpf_pole', hpf_pole, ...
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
tbl_slc     = tbl_impale(neuron_idx, :);
fprintf('--> NOTE: TBL_SLC contains stimuli of %g sec ONLY!\n', unique(tbl_slc.duration_sec));



%% Load PSTH responses
S_list = load.response( tbl_slc, duration_sec );
 
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



%% Save the data
'SAVE the analysis!'
fn.path = load.path_to_data('data');
fn.file = sprintf('data_SU_(%s)_bw(%g)_fbands(%d)_win(%g)ms_spec(%s).mat', ...
    date, binwidth, n_bands, win_size_ms, spectrogram_type);
fn.save = fullfile(fn.path, fn.file);

fprintf('\nSaving data at:\n');
disp(fn)
save(fn.save, 'H', 'tbl_impale', 'neuron_list', 'spec_st', 'S_list');







