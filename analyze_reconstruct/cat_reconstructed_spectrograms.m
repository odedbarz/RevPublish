%
% cat_reconstructed_spectrograms.m
%
% Description:
% Concatenate reconstructed spectrograms into one big matrix.
%
%

clc
% close all
% clear all

fignum = 11;
verbose = 1;

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');


%% Load the measurement's table
% if ~exist('tbl_impale', 'var')
% Save time, load this table only once 
tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');

% Valid list of spiking measurements
%spk_list = find(arrayfun(@(X) str2double(X), tbl_impale.SPK));
spk_list = find( tbl_impale.SPK );
% end
    
drr = get_DRR_list_and_indices;


 
 %% Define the file information
finfo.type = 'MUA'   	% {'SU', 'MUA'}
aux.cprintf('Keywords', '\n-> Data type: *** %s ***\n', finfo.type);

% trainDRR = [3, 5];
% train_DRR_labels = [];
% train_DRR_idx = [];
% for k = 1:length(trainDRR)
%     train_DRR_labels = [train_DRR_labels, drr.labels{trainDRR(k)}, ', '];
%     train_DRR_idx    = [train_DRR_idx, num2str(trainDRR(k)), ', '];
% end
% train_DRR_labels = train_DRR_labels(1:end-2);
% train_DRR_idx = train_DRR_idx(1:end-2);
% aux.cprintf('Keywords', '-> trainDRR: [%s] (%s)\n\n', train_DRR_labels, train_DRR_idx);

if strcmpi('MUA', finfo.type)    
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    %{
    finfo.n_neurons  = [10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '13-Jul-2020';        
    finfo.path       = '../.data/reconstruct/MUA_(13-Jul-2020)_bw(5)ms_fbands(30)_lags(100)ms_cau(0)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_fbands(30)_splits(12)_lags(100)ms_cau(0)_trainDRR(%s).mat';
    %}
        
    % ***
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    % %{
    finfo.n_neurons  = 100 %[1 5 10 25 50 100 150] % 241];  	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '04-Sep-2020';        
    finfo.path       = '../.data/reconstruct/MUA_(04-Sep-2020)_bw(5)ms_algo(svd)_lags(30)ms_cau(0)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_algo(svd)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(%s).mat';
    %}

    % BINWIDTH      : 10 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : multitaper   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '10-Nov-2020';        
    finfo.path       = '../.data/reconstruct/MUA_(10-Nov-2020)_bw(10)ms_fb(50)_lags(30)ms_cau(0)_spec(taper)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(10)ms_fbands(50)_lags(30)ms_cau(0)_trainDRR(%d)_spec(taper_mu5)';
    %}
        
    % BINWIDTH      : 10 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : multitaper   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '10-Nov-2020';        
    finfo.path       = '../.data/reconstruct/MUA_(10-Nov-2020)_bw(10)ms_fb(30)_lags(30)ms_cau(0)_spec(taper)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(10)ms_fbands(30)_lags(30)ms_cau(0)_trainDRR(%d)_spec(taper_mu2)';
    %}
    
    % BINWIDTH      : 10 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : STFT   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '10-Nov-2020';        
    finfo.path       = '../.data/reconstruct/MUA_(10-Nov-2020)_bw(5)ms_fb(50)_lags(30)ms_cau(0)_spec(stft)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_fbands(50)_lags(30)ms_cau(0)_trainDRR(%d)_spec(stft)';
    %}

    % BINWIDTH      : 25 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : gammatone   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '11-Nov-2020';        
    finfo.path       = '../.data/reconstruct/MUA_(11-Nov-2020)_bw(25)ms_fb(30)_lags(30)ms_cau(0)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(25)ms_fbands(30)_lags(30)ms_cau(0)_trainDRR(%d)';
    %}
 
    
    
elseif strcmpi('SU', finfo.type)
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    finfo.n_neurons  = 25 %[10 25 50 103]; % 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '10-Jul-2020';        
    finfo.path       = '../.data/reconstruct/SU_(10-Jul-2020)_bw(5)ms_fbands(30)_lags(100)ms_cau(0)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_fbands(30)_splits(12)_lags(100)ms_cau(0)_trainDRR(%d).mat';
    
else
    error('-> Unrecognized measurement type (%s)!!',  finfo.type);
end
n_neuron_list = length(finfo.n_neurons);



% === Load first batch to get preliminary info
% finfo.fn = sprintf(fn_emplate, finfo.type, finfo.date, finfo.n_neurons(1), finfo.binwidth, ...
%     finfo.n_bands, finfo.win_size_ms, finfo.n_splits, finfo.trainDRR);
finfo.fn = sprintf(fn_emplate, finfo.type, finfo.date, finfo.n_neurons(1), finfo.trainDRR);
warning off
dummy    = load(fullfile(finfo.path, finfo.fn), 'obj_list', 'spec_st', 'splits');
obj_list = dummy.obj_list;
spec_st  = dummy.spec_st;
splits   = dummy.splits;
warning on

% Make sure that the DRY index is the right one used in the loaded file!
%assert(drr.dry == obj_list{1,1}.train_drr);

% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs          = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)
n_bands     = spec_st.n_bands; 
n_time      = spec_st.n_time; 
t           = spec_st.t;            % (sec)


%% Concatenate all the spectrograms and the reconstructions
n_splits = size(obj_list, 2);
drr_to_use = drr.dry;

Sft  = spec_st.Sft{drr_to_use};
Shat = nan(n_bands, n_time);

for k = 1:n_splits
    ix = k == splits.idx;
    Shat(:,ix) = obj_list{drr_to_use, k}.X_est;
end

gof = goodness(Sft, Shat);  

figure(11);
clf;
ax = subplot(2,1,1);
imagesc(t, 1e-3*f, Sft);
set(gca, 'YDir', 'normal', 'XTickLabel', '');
title(sprintf('CC(stim-to-reconstruct): %g', gof.CC));
ylabel('Frequency (kHz)');

ax(2) = subplot(2,1,2);
imagesc(t, 1e-3*f, Shat);
set(gca, 'YDir', 'normal');
xlabel('Time (sec)');
ylabel('Frequency (kHz)');

linkaxes(ax);



%% Save the spectrogram & the reconstructed spectrogram
binwidth    = spec_st.binwidth;     % (ms)
fs          = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)
n_bands     = spec_st.n_bands; 
n_time      = spec_st.n_time; 
t           = spec_st.t;            % (sec)

fn = sprintf('spec_and_reconst_bw(%d)_bands(%d)_units(%d)',...
    binwidth, n_bands, finfo.n_neurons);
save(fn, 'binwidth', 'fs', 'f', 'win_size_ms', 'n_bands', 'n_time', 't', ...
    'Sft', 'Shat');







