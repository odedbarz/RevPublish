function [stim_list, spec_list, meas_list, data_st] = ...
    stimulus_and_spectrogram(T, varargin)
%
% function [T, stim_list, spec_list, meas_list, data_st] = ...
%   load.stimulus_and_spectrogram(neurons, T, varargin)
% 

import load.*


%% Set the inputs
p = inputParser;


% (table) A table that holds all relevant information about the neurons in
%         the database.
addRequired(p, 'T', @istable);            



% (1xn) neuron(s) to load. The type of the loaded stimuli is defined by the
% required neurons. If NEURONS are not defined empty, all available 
% (usually 2, 36 & 40 sec) STIMULI durations will be loaded 
addOptional(p, 'neurons', [], @isnumeric);

% only for spectrogram_type== MULTITAPER
addOptional(p, 'nw', 2, @isnumeric);         

% # of bins along the frequency domain of the spectrogram
addOptional(p, 'n_bands', 30, @isnumeric);    % (1x1)      

% The frequency scale in the spectrogram
addOptional(p, 'f_scale', 'lin', @isstr);       % {'lin', 'log'} 

% The stimulus sampling rate. If the loaded sampling rate is different
%   (e.g., 100kHz), the function will resample the signal
%   * 16kHz was chosen because this is the sampling rate of the WAV files in 
%   the TIMIT corpus.
addOptional(p, 'stim_fs', 16e3, @isstr);       % (1x1. Hz)  

% Spectrogram time-axis sampling rate
% addOptional(p, 'spec_fs_time', 1e3, @isstr);  % (1x1. Hz)  
addOptional(p, 'binwidth', 1, @isnumeric);  % (1x1 ms)  

addOptional(p, 'lowfreq', 250, @isnumeric);     
addOptional(p, 'highfreq', 8000, @isnumeric);     
addOptional(p, 'win_size_ms', 25, @isnumeric);     

% The spectrogram method to use
addOptional(p, 'spectrogram_type', 'matlab', @isstr);    

% If TRUE, performs spectral derivation (difference along the y-axis)
addOptional(p, 'spectral_diff', false, @(x) isnumeric(x) || islogical(x));  

% If TRUE, performs high-pass filtering along the time domain of each 
% frequency channel in the spectrogram.
%
% From Dan Ellis code of find_landmarks.m:
%    ... The high-pass filter applied to the log-magnitude
%        envelope, which is parameterized by the position of the
%        single real pole.  A pole close to +1.0 results in a
%        relatively flat high-pass filter that just removes very
%        slowly varying parts; a pole closer to -1.0 introduces
%        increasingly extreme emphasis of rapid variations, which
%        leads to more peaks initially.
addOptional(p, 'hpf_pole', nan, @(x) isnumeric(x));    


addOptional(p, 'verbose', [], @isnumeric);      % (1x1) write things to the command line?
addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)


parse(p, T, varargin{:});


neurons         = p.Results.neurons;  
nw              = p.Results.nw;  
n_bands         = p.Results.n_bands;
stim_fs         = p.Results.stim_fs;
f_scale         = p.Results.f_scale;
spectrogram_type= p.Results.spectrogram_type;
binwidth        = p.Results.binwidth;      
lowfreq         = p.Results.lowfreq;       
highfreq        = p.Results.highfreq;       
win_size_ms     = p.Results.win_size_ms;   
spectral_diff   = p.Results.spectral_diff;  
hpf_pole        = p.Results.hpf_pole;  

verbose         = p.Results.verbose;       
fignum          = p.Results.fignum;       


if isempty(hpf_pole), hpf_pole = nan; end

%% Init. output parameters
if isempty(neurons)
    % Find all stimulus DURATIONs available in the table T
    durations_all = unique(T.duration_sec);
    assert(~any(isnan(durations_all)), '--> [+load.stimulus_and_spectrogram]: check out your table (excel file)!');
    
    neurons = nan(1, length(durations_all));
    for k = 1:length(durations_all)
        duration_k = durations_all(k);
        neurons(k) = find(T.duration_sec == duration_k, 1, 'first');        
    end    
end
n_neurons = length(neurons);
stim_list = cell(1, n_neurons);
spec_list = cell(1, n_neurons);
meas_list = cell(1, n_neurons);
data_st   = struct();


for kk = 1:n_neurons
    neuron_k = neurons(kk);

    % Set the k'th desired measurement to load
    meas.ID         = T.ID{neuron_k};   
    meas.measType   = T.measType{neuron_k};   
    meas.session    = T.session(neuron_k);    
    meas.unit       = T.unit(neuron_k);       
    meas.measNum    = T.measNum(neuron_k);    
    meas.syncchan   = 0;
    meas.spikechan  = T.spikechan(neuron_k);  
    meas.CF         = T.CF(neuron_k);          % (Hz)
    
    
    %% Load measurement
    [S, n_sync, spikecount, session_fn] = load.session(meas);
    meas.name = sprintf('%s (spikechan: %d)', session_fn, meas.spikechan);

    
    %% Load the k'th STIMULUS 
    clear stim
    suffix       = S.measParam.Suffix;
    duration_sec = S.stimChans{1}.Source.numTokens;
    fn           = get_stimulus_info(suffix, duration_sec);
    if isempty(fn.template)
        return;
    end

    stim.fs = stim_fs;     % (Hz)
    stim = load.stimuli(stim.fs, fn.path, fn.template);

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
    % Calc. spectrograms for all columns\measurements
    Sft = cell(1, size(stim.Y, 2));
    for k = 1:size(stim.Y, 2)
        [Sft{k}, spec_st] = spec.spectrogram(stim.Y(:,k), stim.fs, ...
            'n_bands', n_bands,...
            'lowfreq', lowfreq,...
            'highfreq', highfreq,...
            'overlap_ratio', 0.80,...
            'binwidth', binwidth,...
            'win_size_ms', win_size_ms, ...
            'nw', nw,...                only for spectrogram_type== MULTITAPER
            'f_scale', f_scale,...
            'db_floor', -100, ...  % (dB)
            'duration_ms', duration_ms,...
            'method', spectrogram_type, ...
            'fignum', k+fignum ...
         );
     
        % Perform derivitive along the frequency domain
        if 1 == spectral_diff
            Sft{k} = diff([Sft{k}; Sft{k}(end,:)], [], 1);
        end
     
        % Performs high-pass filtering along the time domain
        if ~isnan(hpf_pole)
            Sft{k} = filter([1 -1], [1 -hpf_pole], Sft{k}')';
        end
    end

    % Add the labels of the spectrograms for later reference
    spec_st.labels = labels(:);
    
    spec_st.Sft = Sft;
    

    %% Save in the lists
    if 1 < n_neurons
        stim_list{kk} = stim;    
        meas_list{kk} = meas;
        spec_list{kk} = spec_st;
    else
        stim_list = stim;    
        meas_list = meas;
        spec_list = spec_st;
    end

end

%% Information about the loaded measurements
valid_table_entries = 1:size(T,1);
[data_st.n_neurons, data_st.n_sessions] = how_many_neurons( T(valid_table_entries,:) );

if verbose
    fprintf('--> # of neurons : %d\n', data_st.n_neurons);
    fprintf('--> # of sessions: %d\n', data_st.n_sessions);   
end

















