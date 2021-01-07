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

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

% Add CochlearModels directory to the path (e.g., Stim2ANF.m).
addpath(genpath('../CochlearModels'));

FigSetup;


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
if ~exist('stim_list', 'var')
    fprintf('--> Load stimuli & spectrograms\n');
    spectrogram_type = 'matlab';      % {['matlab'], 'stft', 'multitaper', 'gammatone'}
    f_scale     = 'log';	% {['lin'], 'log', 'erb'}
    n_bands     = 30       % (1x1) # of bins along the frequency domain of the spectrogram
    binwidth    = 10       % (ms) binwidth of the resulted spectrogram 
    win_size_ms = 50        % (ms) temporal window size over which to calc the spectrogram; 
                            %      'gammatone' filterbanks do not use it!
    lowfreq     = 250;      % (Hz)
    highfreq    = 8000;     % (Hz) %8800;     
    nw          = [];       % applies only for SPECTROGRAM_TYPE = 'multitaper'
    neurons     = 1;        % neuron #1 (#115) is for stimulus of 36 sec (40 sec)
    spectral_diff= 0;       % (logical) perform derivative (DIFF) along the frequency domain
    hpf_pole     = 0.75;
    
    % STIM_LIST & SPEC_LIST contains all stimuli duration (i.e., the 36 s &
    % 40 s stimuli)
    [stim_list, spec_list] = load.stimulus_and_spectrogram(tbl_impale, ...
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
    
    if iscell(spec_list)
        assert(spec_list{1}.binwidth == spec_list{2}.binwidth);
        
        % All available stimulus durations
        duration_all_sec = cellfun(@(X) X.info.Duration , stim_list, 'UniformOutput', 1);
        
        % Plot spectrogram
        %spec.plot_spectrogram(spec_list{1}.t, spec_list{1}.f, spec_list{1}.Sft{3});
    end
        
end

fs = 1/(1e-3*binwidth);         % (Hz)

% ALways load stimuli with the same duration
'Loads stimuli with the SAME duration'
duration_sec= 36;   % (sec) stimulus duration to use 
duration_ms = units.sec2ms( duration_sec );

% Pick the right stimulus parameters for the current loaded measurement
if iscell(spec_list)
    idx_stimuli = find(duration_all_sec == duration_sec);
    stim_st     = stim_list{idx_stimuli};
    spec_st     = spec_list{idx_stimuli};       
else
    stim_st     = stim_list;
    spec_st     = spec_list;       
end

% Plot spectrogram
spec.plot_spectrogram(spec_st.t, spec_st.f, spec_st.Sft{3});


%% Selected neurons table
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
n_drr 	= length(drr.sortby);  
fs_raw  = 10e3;  % (Hz) that's the frequency that I save in the raw table files
fs_mua  = fs;    % (Hz) that's the new downsampled frequency of the MUA data
n_smp   = fs_mua * duration_sec;

neuron_list = find(tbl_impale.duration_sec == duration_sec);
tbl_slc     = tbl_impale(neuron_list, :);
n_rows      = height(tbl_slc);
H           = nan(n_smp, n_drr, n_rows);
H_labels    = zeros(n_rows, n_drr);



%% Prepare the MUA measurements to match the spectrograms
aux.vprint(verbose, '\n-> Starting main loop...\n');
for k = 1:n_rows
    tbl_neuron_k = tbl_slc(k,:);
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
    assert(~isempty(dir(fn.raw)), '--> ERROR: can''t find this file!!\n');
    raw_st = load(fn.raw, '-mat');
    assert(fs_raw == raw_st.sr);

    % Calculate MUA 
    [Hk, ~,  tbl_labels] = calc_raw_means(raw_st.tbl, raw_st.sr, binwidth, duration_sec, tbl_neuron_k.spikechan, 'MUA');    
    assert( size(Hk,1) == duration_sec*1/(1e-3*binwidth), '--> Check out the spectrogram BINWIDTH & the sampling rate of the MUA data!' );
    
    % Save Hk into the big matrix H
    %   The following code sorts the available measurements and save them in 
    %   the right place in H
    cn = 0;     % counter for Hk (some measurements might skip entries, 
                % e.g,  q:1 -> 2 -> 3 -> 5 due to missing
                % measurements).
    for q = 1:height(tbl_labels)
        idx_dist = drr.dist == tbl_labels.Dist(q);
        idx_revb = drr.revb == tbl_labels.Reverb(q);
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


%% Save the data
% %{
'SAVE the analysis!'
fn.save = sprintf('../_data/data_MUA_(%s)_bw(%g)_fbands(%d)_win(%g)ms_spec(%s)_HPF(%g).mat',...
    date, binwidth, n_bands, win_size_ms, spectrogram_type, hpf_pole);
    
save(fn.save, 'H', 'H_labels', 'tbl_impale', 'neuron_list', 'spec_st', 'stim_st');
%}




