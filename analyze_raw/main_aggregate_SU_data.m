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

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

% Add CochlearModels directory to the path (e.g., Stim2ANF.m).
addpath(genpath('../CochlearModels'));

FigSetup


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the measurement's table
if ~exist('tbl_main', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Valid list of spiking measurements
    %spk_list = find(arrayfun(@(X) str2double(X), tbl_impale.SPK));
    %spk_list = find( tbl_impale.SPK );
end
    


%% Load stimuli & measurement data
% % ---------------------------------------------------------
% '*** spk_list ***'
% neurons = spk_list;
% % ---------------------------------------------------------

% if ~exist('stim_list', 'var')
fprintf('--> Load stimuli & spectrograms\n');
spectrogram_type = 'gammatone';      % {['matlab'], 'stft', 'multitaper', 'gammatone'}
f_scale     = 'log';	% {['lin'], 'log', 'erb'}
n_bands     = 30;       % (1x1) # of bins along the frequency domain of the spectrogram
binwidth    = 1;        % (ms) binwidth of the resulted spectrogram 
win_size_ms = nan;      % (ms) temporal window size over which to calc the spectrogram
lowfreq     = 250;      % (Hz)
highfreq    = 8000;     % (Hz)

% STIM_LIST & SPEC_LIST contains all stimuli duration (i.e., the 36 s &
% 40 s stimuli)
[stim_list, spec_list] = load.stimulus_and_spectrogram(tbl_impale, ...
    'spectrogram_type', spectrogram_type, ...
    'binwidth', binwidth,...
    'lowfreq', lowfreq, ...         % (Hz)
    'highfreq', highfreq, ...       % (Hz)
    'win_size_ms', win_size_ms, ...
    'n_bands', n_bands, ...
    'f_scale', f_scale, ... {'lin', 'log'}
    'nw', [], ...          (default: 1.4) only for spectrogram_type == MULTITAPER
    'fignum', [] ...
    );
assert(spec_list{1}.binwidth == spec_list{2}.binwidth);

% All available stimulus durations
duration_all_sec = cellfun(@(X) X.info.Duration , stim_list, 'UniformOutput', 1);

% Plot the spectrogram of the DRY condition
spec.plot_spectrogram(spec_list{1}.t, spec_list{1}.f, spec_list{1}.Sft{3});
% end



fs = 1/(1e-3*binwidth);         % (Hz)

% ALways load stimuli with the same duration
'Loads stimuli with the SAME duration'
duration_sec= 36;   % (sec) stimulus duration to use 
duration_ms = 1e3*duration_sec;     % sec --> ms


% Pick the right stimulus parameters for the current loaded measurement
idx_stimuli = find(duration_all_sec == duration_sec);
stim_st     = stim_list{idx_stimuli};
spec_st     = spec_list{idx_stimuli};  



%% Selected neurons table
% Select the desired stimulus (36 sec or 40 sec)

neuron_idx  = tbl_impale.SPK & (tbl_impale.duration_sec == duration_sec);
neuron_list = find(neuron_idx);
tbl_slc     = tbl_impale(neuron_idx, :);
fprintf('--> NOTE: TBL_SLC contains stimuli of %g sec ONLY!\n', unique(tbl_slc.duration_sec));



%% Load PSTH responses
S_list = load.response( tbl_slc, duration_sec );

% if verbose
%     fprintf('\n--> STIMULUS...\n');
%     fprintf('\t--> Sampling rate (fs): %g kHz\n', 1e-3*stim_st.fs);
%     fprintf('\t--> Duration          : %g (sec)\n', 1e-3*duration_ms);
%     fprintf('\t--> All labels:\n')
% end
 
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
fn.save = sprintf('../.data/data_SU_(%s)_bw(%g)_fbands(%d)_win(%g)ms_spec(%s).mat',...
    date, binwidth, n_bands, win_size_ms, spectrogram_type);
    
save(fn.save, 'H', 'tbl_impale', 'neuron_list', 'spec_st', 'S_list');







