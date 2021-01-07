function [S_list, data_st] = response(tbl_data, duration_sec, varargin)
%
% function [S_list, data_st] = response(tbl_data, duration_sec, varargin)
%
% Inputs:
%   tbl_data
%   duration_sec
%
% Description:
% The function loads a list of Impale measurements. It always loads the same 
% DURATION (sec) at a time.


import load.*


%% Set the inputs
p = inputParser;


% (table) A table that holds all relevant information about the neurons in
%         the database.
addRequired(p, 'tbl_data', @(x) istable(x) || isempty(x));           
addRequired(p, 'duration_sec',  @isnumeric);           

% Duration to load, in seconds
%addOptional(p, 'duration_ms', 36e3, @isnumeric);     % (ms)       

addOptional(p, 'verbose', [], @isnumeric);      % (1x1) write things to the command line?
% addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)

% parse(p, tbl_data, binwidth, varargin{:});
parse(p, tbl_data, duration_sec, varargin{:});

%tbl_data   = p.Results.tbl_data;
%duration_ms= p.Results.duration_ms;
verbose    = p.Results.verbose;       


duration_ms = units.sec2ms(duration_sec);


%% Load DATA into a table
n_neurons = size(tbl_data, 1);
S_list    = cell(1, n_neurons);
data_st   = [];
loaded_entries = false(n_neurons, 1);


%% Init. analysis table
for k = 1:n_neurons    
    % Load the k'th neuron
    meas_k = tbl_data(k, :);
    
    % Make these available in the workspace
    spikechan = meas_k.spikechan;

    % Load measurement
    [S, n_sync, spikecount, session_fn] = load.session(meas_k);
    meas_k.name = sprintf('%s (spikechan: %d)', session_fn, meas_k.spikechan);
          
    % Always load only ONE duration type
    duration_kk_sec = S.stimChans{1}.Source.numTokens;
    if norm(duration_kk_sec-1e-3*duration_ms) > 5*eps
        aux.vprint(verbose, '--> [load.response]: measurement %d was skipped (duration mismatch)!\n', k);
        continue;
    end    
    S.info.duration_ms = duration_ms;
    S.info.spikechan   = spikechan;   
    S.info.n_sync      = n_sync;   
    S.info.spikecount  = spikecount;   
    
    S_list{k} = S;
    loaded_entries(k) = true;
end
 

% Remove all empty entries
idx_nonempty= cellfun(@(X) ~isempty(X), S_list);
S_list      = S_list(idx_nonempty);

% Get back a table with the loaded neurons (remove the others)
tbl_loaded  = tbl_data(loaded_entries, :);

% More information about the loaded measurements
valid_table_entries = 1:size(tbl_loaded,1);
[data_st.n_neurons, data_st.n_sessions] = ...
    how_many_neurons( tbl_loaded(valid_table_entries,:) );
data_st.tbl_loaded = tbl_loaded;

if verbose
    fprintf('--> # of neurons : %d\n', data_st.n_neurons);
    fprintf('--> # of sessions: %d\n', data_st.n_sessions);   
end

















