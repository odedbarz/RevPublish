% plot_CC_RMD_MG.m
%
% Plots one selected unit (single unit or multiunit activity) for the 
% presentation.

clc
fignum = 11;
verbose = 1;
%addpath('../');
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



%% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs          = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
n_drr       = drr.n_drr;
drr_idx     = drr.sortby(1:n_drr);



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
CCs = nan(1, n_drr);    % CCs of responses
CCs_std = nan(1, n_drr);    % CCs of responses
for k = 1:n_drr
    Sk = spec_st.Sft{k};
    
    % CCs: correlation between SRY & DRR stimuli
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



%% Analyze for RMD, MG, CCer, and KURTOSIS
% RMD: response modulation depth
% MG : modulation gain
% CCer: correlation coefficient between DRY envelope (on BF) and resposnes
n_units = size(H_valid,3);

%     '####### DEBUG ######'
%     n_units = 100

% Select data for analysis
H           = squeeze( H_valid(:,:, 1:n_units) );
scores.CCer = nan(n_units, n_drr);       % % CC(response dry & response)
scores.MDr  = nan(n_units, n_drr);
scores.MG   = nan(n_units, n_drr);
scores.skew = nan(n_units, n_drr);
scores.ku   = nan(n_units, n_drr);

% tbl_valid = tbl_impale( data.slc.unit_used, : );
% assert( 0==nnz(tbl_valid.CF==0) && 0 == nnz(isnan(tbl_valid.CF)) );


for n = 1:n_units
    % Option #1:
    % Stimulus envelope; find the closest frequency band to the neuron's CF 
    [~, idx_bf] = min(abs(f - tbl_strf.bf(n)));
    
    % Option #2
    % Stimulus envelope: best correlation between envelope and response
    %'####### DEBUG ######'
    %[~, idx_bf] = max(spec_st.Sft{drr.dry} * H(:,drr.dry,n));
    %idx_bf = 30-idx_bf+1;
    
    % Option #3:
    % Random envelope
    %idx_bf = randi(length(f));
    
    
    % Signal Envelope
    yenv = spec_st.Sft{ drr.dry }(idx_bf,:); 

%     corrcoef(yenv, H(:,1,n))
    
    for rv = 1:n_drr
        % Response MD (modulation depth)
        scores.RMD(n,rv) = MD(H(:,rv,n), 1);
        
        % Envelope's MD 
        % Option #1: using the spectrogram
        % %{
        %yenv = spec_st.Sft{rv}(idx_bf,:); 
        cfloor = spec_st.db_floor;              % The spectrogram (Sft) amplitudes 
        cmax = spec_st.max_Sft;                 % are of log scale, so we need to 
        yenv_10 = cmax*10.^((yenv + cfloor)/20); % revert it
        MD_env = MD( yenv_10 );
        %}
        % Option #2: using the stimulus (more computation time)
        %{
        stim_env_rv = calc_stimulus_envelope(stim_st.Y(:,rv), tbl_valid.CF(n), fs, stim_st.fs);
        MD_env = MD( stim_env_rv, 1 );
        %}
        
        % MG (modulation gain) = RMD/MDenv
        scores.MG(n,rv) = db(scores.RMD(n,rv)./MD_env);
        
        % Kurtosis & Skewness
        scores.skew(n,rv) = skewness( H(:,rv,n) );
        scores.kr(n,rv) = kurtosis( H(:,rv,n) );

        % CCer: DRY(response)-to-DRR(response)
        %dummy = corrcoef(H(:,drr.dry,n), H(:,rv,n));
        dummy = corrcoef(yenv, H(:,rv,n));
        scores.CCer(n,rv) = dummy(1,2);
        
     end
end


% Plot RMD\MG\CC\Kurtosis
%{
figure(30+fignum);
clf;

fontsize = 20;
markersize = 30;

ax = subplot(4,1,1);
plot_dotbox(scores.RMD(:,drr_idx), 'labels', drr.labels(drr_idx));
ylabel(aux.ctitle('RMD', '$(\sqrt{2}\sigma/\mu)$'));
set(gca, 'XTickLabel', '');
xlabel('');
title(ax(1), sprintf('%s', data_type));

ax(2) = subplot(4,1,2);
plot_dotbox(scores.MG(:,drr_idx), 'labels', drr.labels(drr_idx));
ylabel('MG (dB)');
set(gca, 'XTickLabel', '');
xlabel('');

ax(3) = subplot(4,1,3);
plot_dotbox(scores.kr(:,drr_idx), 'labels', drr.labels(drr_idx));
ylabel('Kurtosis');
ylim([min(ylim), 15]);
set(gca, 'XTickLabel', '');
xlabel('');
ylim([1.5, 10]);

ax(4) = subplot(4,1,4);
plot_dotbox(scores.CCer(:,drr_idx), 'labels', drr.labels(drr_idx));
ylabel('$CC$');

hold on
plth = plot(CC_broadband_envelope, 'dk:', 'MarkerSize', 0.4*markersize);
plth(2) = plot(CCs(drr_idx), 'sk:', 'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k');
hold off
legend(plth, {'CCe broadband envelope', 'CCs stimulus'});
xlabel('Direct to Reverberation Ratio (dB)');
set(ax, 'FontSize', fontsize);
aux.abc(ax, 'fontsize', 2*fontsize, 'location', 'northwestoutside');
%}


%%
% CCer: correlation coefficients between bandpass-envelope and response
figure(10+fignum);
% clf;

markersize = 50;
fontsize = 32;

switch upper(data_type)
    case 'SU'
        ax = subplot(1,2,1);        
        plot_dotbox(scores.CCer(:,drr_idx), 'labels', drr.labels(drr_idx));
        %ylabel('$CC_{er}$');
        ylabel('CC (Envelope-to-Response)');
        
    case 'MUA'
        ax = subplot(1,2,2);
        plot_dotbox(scores.CCer(:,drr_idx), 'labels', drr.labels(drr_idx));
        ylabel('CC(Envelope-to-Response)');
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
title(sprintf('%d %ss', n_units, data_type));
set(ax, 'FontSize', fontsize);
ylim([0, 1.2]);

hold on
plth = plot(CC_broadband_envelope, 'dk:', 'MarkerSize', 0.4*markersize);
plth(2) = plot(CCs(drr_idx), 'sk:', 'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k');
hold off
legend(plth, {'CC (Broadband Envelope)', 'CC (Stimulus DRY-to-DRR)'});
xlabel('DRR');
set(ax, 'FontSize', fontsize);
aux.abc(ax, 'fontsize', 2*fontsize, 'location', 'northwestoutside');











