%
% analyze_best_set.m
%
% Description:
% Analyze various "best unit set" scenarios.
%

clc
% close all
% clear all

fignum = 11;
verbose = 1;

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup;


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the measurement's table
if ~exist('tbl_impale', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Valid list of spiking measurements
    %spk_list = find( tbl_impale.SPK );
end
    
drr = get_DRR_list_and_indices;
n_drr = 5;  % N DRRs TO USE



%% Load data
% !!! NOTE !!!
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
% Load:
% 
% 'H', 'tbl_impale', 'spec_st', 'stim_st'
%
data_type    = 'MUA'    % {'SU', MUA'}
fn.load.path = '../.data';
switch upper(data_type)
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        
        %fn.load.file = 'data_SU_(09-Jul-2020)_bw(10)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
        fn.load.file = 'data_SU_(10-Jul-2020)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        %fn.load.file = 'data_SU_(26-Aug-2020)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
        
    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        
        %fn.load.file = 'data_MUA_(09-Jul-2020)_bw(10)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
        fn.load.file = 'data_MUA_(10-Jul-2020)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
        %fn.load.file = 'data_MUA_(23-Jul-2020)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone)';
        %fn.load.file = 'data_MUA_(03-Aug-2020)_bw(5)_fbands(60)_win(NaN)ms_spec(gammatone)';
        %fn.load.file = 'data_MUA_(03-Aug-2020)_bw(5)_fbands(30)_win(10)ms_spec(stft)';
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data    = load(fn.load.fullfile);
spec_st = data.spec_st;


aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);
aux.vprint(verbose, '-> data_type: %s\n', data_type);



%% Get the valid measurements\columns
duration_sec = 36;  % (sec) stimulus duration to use

% Select all measurement with the desired duration
slc.duration = tbl_impale.duration_sec == duration_sec;

% Select all measurements with that contains all (5) sessions
% slc_valid_neurons = ~arrayfun(@(I) 5 > nnz(~isnan(sum(H(:,:,I))) | 0 == sum(H(:,:,I))), 1:size(H,3));
slc.valid_neurons = 5 <= sum( ~isnan( squeeze(sum(data.H,1)) ), 1)';
slc.valid_neurons = slc.valid_neurons(:);

% Indices of both boolean conditions
slc.valid_neuron_idx = slc.valid_neurons & slc.duration(data.neuron_list);


% # of units to use in this script
units = 1:nnz(slc.valid_neuron_idx); 

% Choosing # of units
H_valid = data.H(:,:,slc.valid_neuron_idx);
H_unit = squeeze( H_valid(:,:,units) );


% A list of all "valid" units to use
slc.neuron_list = data.neuron_list(slc.valid_neuron_idx);
slc.neuron_list = slc.neuron_list(units);  % truncate the list



%%
iscausal    = 0;                    % use reconstructed causal filters? 
lags_ms     = 30                    % (ms) maximum lags
binwidth    = spec_st.binwidth;     % (ms)
n_bands     = spec_st.n_bands;
win_size_ms = spec_st.win_size_ms; 
% algo_type   = 'regression';     % {'regression', 'asd', 'svd'}
algo_type   = 'svd';     % {'regression', 'asd', 'svd'}
jk_flag     = 1;

% Computation parameters (Reconstuction object)
gpu_flag    = 0;    %0 < gpuDeviceCount;   % use the GPU, if available

if verbose
    aux.cprintf('UnterminatedStrings', '\n    Data:\n');
    aux.cprintf('UnterminatedStrings', '--> data_type   : %s\n', data_type);
    aux.cprintf('UnterminatedStrings', '--> duration_sec: %g ms\n', duration_sec);
    aux.cprintf('UnterminatedStrings', '    Reconstruction:\n');
    aux.cprintf('UnterminatedStrings', '--> causality   : %d\n', iscausal);
    aux.cprintf('UnterminatedStrings', '--> lags_ms     : %g ms\n', lags_ms);
    aux.cprintf('UnterminatedStrings', '--> binwidth    : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands     : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms : %g ms\n', win_size_ms);
    aux.cprintf('UnterminatedStrings', '--> algo_type   : %s ms\n', algo_type);
    aux.cprintf('UnterminatedStrings', '--> is_jackknife: %d\n', jk_flag);
end




%% Split the WAV file into chunks of the various speakers
% Splits the overall stimulus into chunks according to the speakers
[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(binwidth, 1e-3*spec_st.duration_ms);

% Make sure that the split indices have a valid length
assert(spec_st.n_time == split_time_idx(end));

% Get the longest time interval
%n_smp_split = max( diff(split_time_idx, [], 2) );






%% *** SECTION #1 ***
%  ******************
% Analyze a given chunk\speaker
%analyze_best.analyze_chosen_chunk;



%% *** SECTION #2 ***
%  ******************
% Check reconstructions as a function of # of units
%analyze_best.reconstruct_over_units;



%% *** SECTION #3 ***
%  ******************
% Check reconstructions as a function of LAGs
analyze_best.reconstruct_over_lags;



%% Save
%{
fn.save.path = '../.data/reconstruct/';
fn.save.file = sprintf('CCvsUnits_%s_(%s)_bw(%g)ms_fbands(%d)_sortby(%s)_lags(%g)ms_cau(%d)_trainDRR(%d)_testDRR(%d)',...
    data_type, date, binwidth, n_bands, sort_algo, lags_ms, iscausal, train_drr, test_drr);
fn.save.fullfile = fullfile( fn.save.path, fn.save.file );

% Save the results for that 
save(fn.save.fullfile, 'gof', 'X_test', 'y_test', 'obj'); 
%}










