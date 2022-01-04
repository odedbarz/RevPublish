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

setup_environment;


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');



%% Load the measurement's table
if ~exist('tbl_main', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Row is a duplicate variable name used by MATLAB
    tbl_impale.Properties.VariableNames{'Row'} = 'neuron';
end
    



%% Load stimuli & measurement data
fprintf('--> Load stimuli & spectrograms\n');
spectrogram_type = 'gammatone';      % {['matlab'], 'stft', 'multitaper', 'gammatone'}
f_scale     = 'erb';	% {['lin'], 'log', 'erb'}
n_bands     = 30      % (1x1) # of bins along the frequency domain of the spectrogram
binwidth    = 5         % (ms) binwidth of the resulted spectrogram 
win_size_ms = nan       % (ms) temporal window size over which to calc the spectrogram; 
                        %      'gammatone' filterbanks do not use it!
lowfreq     = 100;      % (Hz)
highfreq    = 8900;     % (Hz) 8900 Hz so that to get as high as possible to near 8k Hz with the 
nw          = [];       % applies only for SPECTROGRAM_TYPE = 'multitaper'
duration_to_load = 36;   % neuron #1 (#115) is for stimulus of 36 sec (40 sec)
spectral_diff= 0;       % (logical) perform derivative (DIFF) along the frequency domain
hpf_pole    = nan;

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

% 
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
    aux.cprintf('UnterminatedStrings', '--> lowfreq         : %g Hz\n', lowfreq);
    aux.cprintf('UnterminatedStrings', '--> highfreq        : %g Hz\n', highfreq);
    aux.cprintf('UnterminatedStrings', '--> spectral_diff   : %d\n', spectral_diff);
end




%% *** Get the MUA data into one big matrix ***
% A 3D matrix to hold all the MUA responses
drr     = get_DRR_list_and_indices;
n_drr 	= drr.n_drr;  
fs_raw  = 10e3;  % (Hz) that's the frequency that I save in the raw table files
fs_mua  = fs;    % (Hz) that's the new downsampled frequency of the MUA data
n_smp   = fs_mua * duration_sec;

% Get all rows (measurements) with desired duration
neuron_list = find(tbl_impale.duration_sec == duration_sec);

tbl_MUA     = tbl_impale(neuron_list, :);
n_rows      = height(tbl_MUA);
H           = nan(n_smp, n_drr, n_rows);
H_labels    = zeros(n_rows, n_drr);

% Add another column for all recorded data
tbl_MUA_all= tbl_MUA;
all_meas    = cell(height(tbl_MUA_all), 1);
tbl_MUA_all= [tbl_MUA_all, table(all_meas)];



%% Prepare the MUA measurements to match the spectrograms
aux.vprint(verbose, '\n-> Starting main loop...\n');
for k = 1:n_rows
    tbl_neuron_k = tbl_MUA_all(k,:);
    S = load.response( tbl_neuron_k, duration_sec );
    if isempty(S)
        aux.vprint(verbose, '--> Skipping measurement (%d/%d)...\n', k, n_rows);
        continue;
    elseif 0==mod(k-1,10)
        aux.vprint(verbose, '--> Loading Imaple structure (%d/%d)\n', k, n_rows);
    end
    S  = S{1};     % cell --> struct

    % Loads the raw data
    fn.raw = fullfile(path_root_raw, S.info.rawDataDir, S.info.rawFileName);
    assert(~isempty(dir(fn.raw)), '--> ERROR: can''t find this file: \n\t %s !!\n', fn.raw);
    raw_st = load(fn.raw, '-mat');
    assert(fs_raw == raw_st.sr);

    % Calculate MUA 
    [Hk, ~,  tbl_kth_meas] = calc_raw_means(raw_st.tbl, raw_st.sr, binwidth, duration_sec, tbl_neuron_k.spikechan, 'MUA');    
    assert( size(Hk,1) == duration_sec*1/(1e-3*binwidth), '--> Check out the spectrogram BINWIDTH & the sampling rate of the MUA data!' );
    
    % Save all single measurements
    tbl_MUA_all.all_meas{k} = tbl_kth_meas;
    
    % Save Hk into the big matrix H
    %   The following code sorts the available measurements and save them in 
    %   the right place in H
    cn = 0;     % counter for Hk (some measurements might skip entries, 
                % e.g,  q:1 -> 2 -> 3 -> 5 due to missing
                % measurements).
    for q = 1:height(tbl_kth_meas)
        idx_dist = drr.dist == tbl_kth_meas.Dist(q);
        idx_revb = drr.revb == tbl_kth_meas.Reverb(q);
        idx_q = find( idx_dist & idx_revb );

        % If there is no such label, move on
        if isempty(idx_q)
            %warning('---> [for q = ...]: can''t find this label!');
            continue;
        end

        cn = cn + 1;

        % Save the loaded measurement in the right place 
        H(:,idx_q,k) = Hk(:,cn);

        % Add a flag of the measurement's type
        H_labels(k, idx_q) = 1;

    end

end


%% Get all rows (measurements) that have full session (all DRR conditions)
valid_neuron_idx = n_drr <= sum(H_labels(:,1:n_drr),2);
tbl_MUA = tbl_MUA(valid_neuron_idx, :);
H = H(:,:,valid_neuron_idx);



%% Save the data
% %{
fprintf('SAVE the analysis!');

fn.path = load.path_to_data('data');
fn.file = sprintf('data_MUA_(%s)_bw(%g)_fbands(%d)_win(%g)ms_spec(%s)', ...
    date, binwidth, n_bands, win_size_ms, spectrogram_type);
fn.save = fullfile(fn.path, fn.file);

fprintf('\nSaving data at:\n');
disp(fn)
save(fn.save, 'H', 'tbl_MUA', 'spec_st', 'stim_st');
save([fn.save, '_ALL-MEAS'], 'tbl_MUA_all', 'spec_st', 'stim_st', '-v7.3');
%}




