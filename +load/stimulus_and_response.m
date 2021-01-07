function [tbl_data, spec_list, psth_list, stim_list, meas_list, data_st] = ...
    stimulus_and_response(neurons, varargin)
%
%   [tbl_data, spec_list, psth_list, stim_list, meas_list] = load.stimulus_and_response(neurons, varargin)
%

import load.*


%% Set the inputs
p = inputParser;

% (1xn) neuron(s) to load. If 0 == neurons, the function outputs just the
%       tbl_data.
addRequired(p, 'neurons', @isnumeric);            


% (table) A table that holds all relevant information about the neurons in
%         the database.
addOptional(p, 'tbl_data', [], @(x) istable(x) || isempty(x));           

% The bin width for the spectrogram and the PSTH(s)
addOptional(p, 'win_size_ms', 2.0, @isnumeric);    % (ms)       


% The rabbit's ID to load
addOptional(p, 'rabbit_ID', 'C74', @isstr);           

% The stimulus sampling rate. If the loaded sampling rate is different
%   (e.g., 100kHz), the function will resample the signal
%   * 16kHz was chosen because this is the sampling rate of the WAV files in 
%   the TIMIT corpus.
addOptional(p, 'stim_fs', 16e3, @isstr);       % (1x1. Hz)  

addOptional(p, 'verbose', [], @isnumeric);      % (1x1) write things to the command line?
addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)


parse(p, neurons, varargin{:});

tbl_data    = p.Results.tbl_data;
win_size_ms = p.Results.win_size_ms;       % (ms) PSTHs bin-width
rabbit_ID   = p.Results.rabbit_ID;
stim_fs     = p.Results.stim_fs;
verbose     = p.Results.verbose;       
% fignum   = p.Results.fignum;       



%% Load DATA into a table
data_st.path	= path_to_data; 
data_st.fn      = 'LabNoteBook.xlsx';
data_st.operator= 'OBZ';
data_st.ID      = rabbit_ID;
data_st.measType= 'Spch';   % {'Spch', 'ns_Spch_fc1kHz', 'ns_Spch_fc4kHz'}
data_st.session = [];     
data_st.unit    = [];      
data_st.measNum = [];     

if isempty(tbl_data)
    tbl_data = load.data_table(data_st);
end

n_neurons = length(neurons);

psth_list = cell(1, n_neurons);
stim_list = cell(1, n_neurons);
spec_list = cell(1, n_neurons);
meas_list = cell(1, n_neurons);


if (0 == n_neurons) || (0 == neurons)
    return;
end


%% Init. analysis table
% n_all_neurons   = size(tbl_data, 1);
% binwidth = 1;      % (ms) PSTHs bin-width
counter  = 0; 


for kk = 1:n_neurons
    neuron_k = neurons(kk);
    
    % Set the k'th desired measurement to load
    meas.measType   = tbl_data.measType{neuron_k};   
    meas.session    = tbl_data.session(neuron_k);    
    meas.unit       = tbl_data.unit(neuron_k);       
    meas.measNum    = tbl_data.measNum(neuron_k);    
    meas.syncchan   = 0;
    meas.spikechan  = tbl_data.spikechan(neuron_k);  
    meas.CF         = tbl_data.CF(neuron_k);          % (Hz)

    % Make these available in the workspace
    spikechan       = meas.spikechan;
    syncchan        = meas.syncchan;

    % Load measurement
    [S, n_sync, spikecount, session_fn] = load_session(meas);
    meas.name = sprintf('%s (spikechan: %d)', session_fn, meas.spikechan);
          
      
        
    %% Load the k'th STIMULUS 
    clear stim
    suffix       = S.measParam.Suffix;
    duration_sec = S.stimChans{1}.Source.numTokens;
    fn           = get_stimulus_info(suffix, duration_sec);
    if isempty(fn.template)
        continue;
    end
        
    stim.fs = stim_fs;     % (Hz)
    stim = load_stimuli(stim.fs, fn.path, fn.template);
    
    % keep record of the WAV file
    stim.fn_template = fn.template;   
    stim.fn_path = fn.path;   

    % Duration of the stimulus, in mili-seconds
    duration_ms = units.sec2ms( stim.info.Duration );
    stim.duration_ms = duration_ms;
    
    % Measurement's labels
    labels = stim.labels;
    stim.labels = labels;
    
    % (cell --> matrix)
    stim.Y = [stim.Y{:}];
    
    if verbose
        fprintf('\n');    
        fprintf('--> (%d) Analyzing SESSION: %s\n', kk, meas.name);
        disp( table(n_sync, spikecount) );
        
        fprintf('\n--> STIMULUS: %g\n', stim.fs);
        fprintf('    ---------\n');
        fprintf('--> * Sampling rate (fs): %g\n', stim.fs);
        fprintf('--> * Duration: %g\n', stim.info.Duration);
        fprintf('--> * Loaded indices\n');
        fprintf('--> * Measurement labels:\n')
        disp(labels);
    end
        
    
    %% ** Spectrogram ** 
    spec.n_bands    = 30;       % default: 30
    spec.lowfreq    = 250;      % (Hz)
    spec.highfreq   = 8000;     % (Hz)
    spec.overlap_ratio = 0.50;  % NOVERLAP ratio (out of the window size) in MATLAB's SPECTROGRAM
    spec.fs_stim    = stim.fs;  % (Hz)
    spec.win_size_ms= win_size_ms;      % (ms) default: 40 ms in Mesgarani et al., 2014 
    spec.f_scale    = 'log';    % {'lin', 'log'}
    spec.labels     = labels(:);
    %spec.fignum     = [];

    % Calc. spectrograms for all columns\measurements
    [spec_st, Sft_all] = arrayfun(@(QQ) sgram(stim.Y(:,QQ), stim.fs, spec), ...
       1:size(stim.Y, 2), 'UniformOutput', false);
    
    spec = spec_st{1};
    spec.Sft = Sft_all;

    
    
    %% ** PSTH ** 
    psth_k = PSTH(S,...
        duration_ms,...
        spec.binwidth,...   % !! set the PSTH bins to equal that of the spectrogram
        spikechan,...
        syncchan );
    
   
    %% Save in the lists
    counter = counter +1;   
    if 1 < n_neurons
        stim_list{counter} = stim;    
        meas_list{counter} = meas;
        psth_list{counter} = psth_k;
        spec_list{counter} = spec;
    else
        stim_list = stim;    
        meas_list = meas;
        psth_list = psth_k;
        spec_list = spec;
    end
    
end

tbl_data = tbl_data(neurons, :);


%% Information about the loaded measurements
valid_table_entries = 1:size(tbl_data,1);
[data_st.n_neurons, data_st.n_sessions] = how_many_neurons( tbl_data(valid_table_entries,:) );

if verbose
    fprintf('--> # of neurons : %d\n', data_st.n_neurons);
    fprintf('--> # of sessions: %d\n', data_st.n_sessions);   
end

















