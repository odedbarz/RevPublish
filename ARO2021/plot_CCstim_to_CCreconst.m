% plot_CCstim_to_CCreconst.m
%
% Plots one selected unit (single unit or multiunit activity) for the 
% presentation.

clc
fignum = 11;
verbose = 1;

addpath('../');

FigSetup;



%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the mea
if ~exist('tbl_impale', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Valid list of spiking measurements
    %spk_list = find( tbl_impale.SPK );
end
    
drr = get_DRR_list_and_indices;
%n_drr = 5;  % N DRRs TO USE


%% META-DATA 
dummy = load('..\.data\stimulus\metadata_(36)_wav_(30-Jun-2020).mat');
tbl_metadata = dummy.tbl_metadata;



%% Load data
% %{
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
        analyze.path = '../.data/Analyze/';
        analyze.filename = 'Analysis_SU_(02-Sep-2020)_units(103)_bw(5)_fbands(30).mat  ';
        load([analyze.path, analyze.filename]);
        
    case 'MUA'
        analyze.path = '../.data/Analyze/';
        analyze.filename = 'Analysis_MUA_(02-Sep-2020)_units(100)_bw(5)_fbands(30).mat';
        load([analyze.path, analyze.filename]);
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
% fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
% data    = load(fn.load.fullfile);
% spec_st = data.spec_st;
% aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);
% aux.vprint(verbose, '-> data_type: %s\n', data_type);
%}



%% Get the valid measurements\columns
%{
duration_sec = 36;  % (sec) stimulus duration to use

% Select all measurement with the desired duration
neuron_duration_list = tbl_impale.duration_sec == duration_sec;

% Select all measurements with that contains all (5) sessions
slc.valid_neuron_idx = 5 <= sum( ~isnan( squeeze(sum(data.H,1)) ), 1)';
slc.valid_neuron_idx = slc.valid_neuron_idx(:);

% Indices of both boolean conditions
slc.valid_neuron_idx = slc.valid_neuron_idx & neuron_duration_list(data.neuron_list);

% A list of all "valid" units to use
slc.valid_neurons = data.neuron_list(slc.valid_neuron_idx);
%n_units = nnz(slc.valid_neuron_idx);

% Choosing # of units
H_valid = data.H(:,:,slc.valid_neuron_idx);
%}


%% Sampling frequency along the time axis
% %{
binwidth    = spec_st.binwidth;     % (ms)
fs          = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
n_drr       = drr.n_drr;
drr_idx     = drr.sortby(1:n_drr);
%}


%% BROAD-BAND CCs
% Load the stimuli structures
dummy          = load(fullfile('../.data/stimulus', 'data_stimuli_duration(36_and_40)sec.mat'));
stim_list      = dummy.stim_list;

% All available stimulus durations
duration_all_sec = cellfun(@(X) X.info.Duration , stim_list, 'UniformOutput', 1);

% Select the desired stimulus DURATION
duration_sec= 36;  % (sec)
duration_ms = 1e3*( duration_sec );
idx_stimuli = find(duration_all_sec == duration_sec);
stim_st     = stim_list{idx_stimuli};

% Calculate the broadband CCs
win_env_lpf_ms = 60;
fs_dwn         = fs;
[R, args]      = broadband_envelopes_CCs(stim_st.Y, stim_st.fs, win_env_lpf_ms, fs_dwn);
CC_broadband_envelope = R(drr.dry, drr.sortby(1:n_drr));



%% Correlation Coefficients between STIMULI
Sdry = spec_st.Sft{drr.dry};

CCs = nan(1, n_drr);    % CCs stimulus DRY-to-DRR
CCs_std = nan(1, n_drr);    % CCs of responses
for k = 1:n_drr
    Sk = spec_st.Sft{k};
    
    % CCs stimulus DRY-to-DRR
    cs_k = diag( zca(Sdry')' * zca(Sk') )/size(Sk,2);
    CCs(k) = mean(cs_k);
    CCs_std(k) = sqrt(mean((cs_k/size(Sk,2) - CCs(k)).^2));    
end



%% Load the STRFs data for the STRF-BEST FREQUENCIES 
% %{
% whos
%   Name               Size                  Bytes  Class     Attributes
% 
%   H_valid         7200x6x103            35596800  double              
%   scores             1x1                  149976  struct              
%   slc                1x1                    1326  struct              
%   spec_st            1x1                15039301  struct              
%   splits             1x1                   58592  struct              
%   tbl_impale       437x20                 640191  table               
%   tbl_strf         103x5                  373423  table               
clear tbl_strf
fn.path = '..\.data\STRFs\';
switch upper(data_type)
    case 'SU'
        fn.name = 'STRF_SU_(25-Aug-2020)_units(103)_bw(5)ms_algo(regression)_fbands(30)_lags(30)ms_cau(1)_trainDRR(3).mat';
        % fn.name = 'STRF_SU_(27-Aug-2020)_units(103)_bw(1)ms_algo(regression)_fbands(30)_lags(20)ms_cau(1)_trainDRR(3)';
        
    case 'MUA'
        fn.name = 'STRF_MUA_(02-Sep-2020)_units(241)_bw(5)ms_algo(regression)_fbands(30)_lags(20)ms_cau(1)_trainDRR(3).mat';
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
        
fn.file2load = fullfile(fn.path, fn.name);
data_strf = load(fn.file2load);
tbl_strf = data_strf.tbl_strf;
assert(data_strf.spec_st.binwidth == data_strf.spec_st.binwidth);
%}



%% Analyze for CCer
n_units = size(H,3);
CCer    = nan(n_units, n_drr); 
for n = 1:n_units
    % Option #1:
    % Stimulus envelope; find the closest frequency band to the neuron's CF 
    [~, idx_bf] = min(abs(f - tbl_strf.bf(n)));
    
    % Option #2
    % Stimulus envelope: best correlation between envelope and response
    %'####### DEBUG ######'
    %[~, idx_bf] = max(spec_st.Sft{drr.dry} * H(:,drr.dry,n));
    %idx_bf = 30-idx_bf+1;
    
    
    % Signal Envelope
    yenv = spec_st.Sft{ drr.dry }(idx_bf,:); 
   
    for rv = 1:n_drr
        % Response MD (modulation depth)
        scores.RMD(n,rv) = MD(H(:,rv,n), 1);
        
        % CCer: DRY(response)-to-DRR(response)
        %dummy = corrcoef(H(:,drr.dry,n), H(:,rv,n));
        dummy = corrcoef(yenv, H(:,rv,n));
        CCer(n,rv) = dummy(1,2);
        
     end
end


%% Plot
fontsize = 28;
markersize = 40;


% CC(DRY, DRR)
plth = plot(1:n_drr, CCs(drr_idx), 'sk:', 'MarkerSize', 0.4*markersize);
arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));

hold on
% CC of broadband envelopes (DRY to DRR)
plth(2) = plot(CC_broadband_envelope, 'dk:', 'MarkerSize', 0.4*markersize);

% Plot CCer(DRY(response)-to-DRR(response)) vs. CC(Sdry-to-\hat{S}drr)
M = [mean(CCer(:,drr.ordered))', scores.mu.CC(drr.ordered,end)];
errbar = [std(CCer(:,drr.ordered))', scores.SE.CC(drr.ordered,end)];
h = bar(M); 
set(gca, 'XTick', 1:drr.n_drr, 'XTickLabel', drr.labels(drr.ordered));
h(1).FaceColor = aux.rpalette(sprintf('new%02d',1));
h(2).FaceColor = aux.rpalette(sprintf('new%02d',2));
hold off












