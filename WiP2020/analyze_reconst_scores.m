%
% analyze_reconst_scores.m
%
% Description:
% This code analyze the reconstructed spectrograms and compares these to
% the stimulus spectrograms (DRY and\or other DRR conditions).
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
%path_root_raw = load.path_to_data('raw');


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
    finfo.path       = '../_data/reconstruct/MUA_(13-Jul-2020)_bw(5)ms_fbands(30)_lags(100)ms_cau(0)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_fbands(30)_splits(12)_lags(100)ms_cau(0)_trainDRR(%s).mat';
    %}
    
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    %{
    finfo.n_neurons  = [10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '04-Sep-2020';        
    finfo.path       = '../_data/reconstruct/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_algo(svd)_fbands(30)_splits(12)_lags(20)ms_cau(0)_trainDRR(%d).mat';
    %}
    
    % *** DEFAULT ***
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    % %{
    finfo.n_neurons  = 100 %[1 5 10 25 50 100 150] % 241];  	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '04-Sep-2020';        
    finfo.path       = '../_data/reconstruct/MUA_(04-Sep-2020)_bw(5)ms_algo(svd)_lags(30)ms_cau(0)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_algo(svd)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(%s).mat';
    %}

    % BINWIDTH      : 10 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : multitaper   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '10-Nov-2020';        
    finfo.path       = '../_data/reconstruct/MUA_(10-Nov-2020)_bw(10)ms_fb(50)_lags(30)ms_cau(0)_spec(taper)/';
    %fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(10)ms_fbands(50)_lags(30)ms_cau(0)_trainDRR(%s)_spec(taper_mu5)';
    fn_emplate       = 'rc_%s_(%s)_units(%d)_bw(10)ms_fb(50)_lags(30)ms_cau(0)_trDRR(%s)_spec(taper_mu5)';
    %}
        
    % BINWIDTH      : 10 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : multitaper   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '10-Nov-2020';        
    finfo.path       = '../_data/reconstruct/MUA_(10-Nov-2020)_bw(10)ms_fb(30)_lags(30)ms_cau(0)_spec(taper)/';
    %fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(10)ms_fbands(30)_lags(30)ms_cau(0)_trainDRR(%d)_spec(taper_mu2)';
    fn_emplate       = 'rc_%s_(%s)_units(%d)_bw(10)ms_fb(30)_lags(30)ms_cau(0)_trDRR(%s)_spec(taper_mu2)';
    %}
    
    % BINWIDTH      : 10 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : STFT   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '10-Nov-2020';        
    finfo.path       = '../_data/reconstruct/MUA_(10-Nov-2020)_bw(5)ms_fb(50)_lags(30)ms_cau(0)_spec(stft)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_fbands(50)_lags(30)ms_cau(0)_trainDRR(%d)_spec(stft)';
    %}

    % BINWIDTH      : 25 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : gammatone   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '11-Nov-2020';        
    finfo.path       = '../_data/reconstruct/MUA_(11-Nov-2020)_bw(25)ms_fb(30)_lags(30)ms_cau(0)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(25)ms_fbands(30)_lags(30)ms_cau(0)_trainDRR(%d)';
    %}
    
    % !! TESTING !!
    % BINWIDTH      : 25 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : gammatone   
    %{
    finfo.n_neurons  = 100; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '11-Nov-2020';        
    finfo.path       = '../_data/reconstruct/MUA_(11-Nov-2020)_bw(5)ms_fb(30)_lags(30)ms_cau(0)/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_fbands(30)_lags(30)ms_cau(0)_trainDRR(%d)';
    %}
    
    % !! MORE than one DRR condition !!
    % BINWIDTH      : 5 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : gammatone   
    %{
    finfo.n_neurons  = 50; %[10 25 50 100 150];  	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '11-Nov-2020';        
    finfo.path       = '../_data/reconstruct/MUA_(11-Nov-2020)_bw(5)ms_fb(30)_lags(30)ms_cau(0)_trDRR(2)/';
    fn_emplate       = 'rct_%s_(%s)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(%s)';
    %}

    % '*** MANIPULATING THE SPECTROGRAM ***'
    %{
    % BINWIDTH      : 5 ms
    % CAUSALITY     : false
    % SPECTROGRAM   : gammatone   
    % %{
    finfo.n_neurons  = 241;  	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '08-Dec-2020';        
    finfo.path       = '../_data/reconstruct/';
    fn_emplate       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(%s)';
    %}

  
    
elseif strcmpi('SU', finfo.type)
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    finfo.n_neurons  = 25 %[10 25 50 103]; % 150];  	% # of units to load from existing files
    finfo.trainDRR   = trainDRR;
    finfo.date       = '10-Jul-2020';        
    finfo.path       = '../_data/reconstruct/SU_(10-Jul-2020)_bw(5)ms_fbands(30)_lags(100)ms_cau(0)/';
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
dummy    = load(fullfile(finfo.path, finfo.fn), 'obj_list');
obj_list = dummy.obj_list;
dummy    = load(fullfile(finfo.path, finfo.fn), 'spec_st');
spec_st  = dummy.spec_st;
warning on

% Make sure that the DRY index is the right one used in the loaded file!
%assert(drr.dry == obj_list{1,1}.train_drr);

% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs          = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)

% OBJ_LIST
% * Dim1: # of DRR cases checked 
% * Dim2: # of  simulations 
[n_drr, n_JK] = size(obj_list);
assert(5 == n_drr);
n_splits = n_JK;
n_bands  = obj_list{1}.n_bands;


if verbose
    aux.cprintf('UnterminatedStrings', '--> binwidth   : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands    : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms: %g ms\n', win_size_ms);
    aux.cprintf('UnterminatedStrings', '--> n_splits   : %d\n', n_splits);
end


%% BROAD-BAND CCs
% Load the stimuli structures
dummy          = load(fullfile('../_data/stimulus', 'data_stimuli_duration(36_and_40)sec.mat'));
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



%% META-DATA 
dummy = load('../_data/stimulus/metadata_(36)_wav_(30-Jun-2020).mat');
tbl_metadata = dummy.tbl_metadata;

tbl_dummy = [array2table(tbl_metadata.fn, 'VariableNames', {'fn'}), array2table([(1:12)', [tbl_metadata.t0, tbl_metadata.t1]/16e3])];
tbl_dummy.Properties.VariableNames(2:4) = {'N', 't0_sec', 't1_sec'};
disp( tbl_dummy );



%% SCORES & STI structures
clear scores scores2 sti sti2
scores.n_drr        = n_drr;
scores.n_JK         = n_JK;
scores.n_neuron_list= n_neuron_list;
scores.n_bands      = n_bands;
scores.n_units      = nan(1, n_neuron_list);
scores.CC           = nan(n_drr, n_JK, n_neuron_list);
scores.mse          = nan(n_drr, n_JK, n_neuron_list);
scores.nmse         = nan(n_drr, n_JK, n_neuron_list);
scores.CCf          = nan(n_drr, n_JK, n_neuron_list, n_bands);

sti.n_drr        = n_drr;
sti.n_JK         = n_JK;
sti.units        = finfo.n_neurons;
sti.n_neuron_list= n_neuron_list;
sti.test         = nan(n_drr, n_JK);        % X_test is always DRY condition!
sti.test_probe   = nan(n_drr, n_JK);        % X_test is always DRY condition!
sti.est          = nan(n_drr, n_JK, n_neuron_list);
sti.est_probe    = nan(n_drr, n_JK, n_neuron_list);

% For comparing between stimuli of various DRR conditions and the
% reconstructed spectrograms
%{
sti2.n_drr        = n_drr;
sti2.n_JK         = n_JK;
sti2.n_neuron_list= n_neuron_list;
sti2.est          = nan(n_drr, n_JK, n_neuron_list);
%}
scores2           = scores;



%%
for n = 1:n_neuron_list
    % Load the data
    %          H: [3600󬝭50 double]
    %      X_est: [30�180󬊂0 double]
    %     X_test: [30�180 double]
    %   obj_list: {5�20 cell}
    %    spec_st: [1�1 struct]
    % tbl_impale: [437�19 table]
    %finfo.fn = sprintf(fn_emplate, finfo.type, finfo.date, finfo.n_neurons(n), finfo.binwidth, ...
    %    finfo.n_bands, finfo.win_size_ms, finfo.n_splits, finfo.trainDRR);
    finfo.fn = sprintf(fn_emplate, finfo.type, finfo.date, finfo.n_neurons(n), finfo.trainDRR);
    aux.cprintf('g', '\n-> Loading <%s>...\n', finfo.fn);
    warning off
    data    = load(fullfile(finfo.path, finfo.fn));
    warning on
    obj_list= data.obj_list;
    spec_st = data.spec_st;
    H       = data.H_units;
    splits  = data.splits;
    
    % # of neurons
    n_neurons = obj_list{1}.n_neurons;

    % Make sure that the loaded data and the requested data have the same
    % parameters
    assert(binwidth    == spec_st.binwidth);
    assert(n_bands     == spec_st.n_bands);
    if ~isempty(spec_st.win_size_ms)
        assert(isnan(win_size_ms)&&isnan(spec_st.win_size_ms) || (win_size_ms == spec_st.win_size_ms));
    end
    assert(n_JK        == splits.n_grps );

    if verbose
        fprintf('\n');
        fprintf('--> trainDRR   : %s\n', finfo.trainDRR);
        fprintf('--> n_neurons  : %d\n', finfo.n_neurons(n));
        fprintf('--> binwidth   : %d\n', binwidth);
        fprintf('--> n_bands    : %d\n', n_bands);
        fprintf('--> win_size_ms: %d\n', win_size_ms);
        fprintf('--> n_drr      : %d\n', n_drr);
        fprintf('--> n_JK       : %d\n', n_JK);
    end
    

    %%
    scores.n_units(n) = n_neurons;
    for jk = 1:n_JK 
        % Cut out the testing chunk of the spectrogram
        % * n_JK == test_grp_number         
        X_dry = spec_st.Sft{drr.dry}(:, jk == splits.idx);    
        
        for rv = 1:n_drr               
            % Cut out the testing chunk of the spectrogram
            % * n_JK == test_grp_number 
            X_est                 = obj_list{rv,jk}.X_est;
            gof                   = goodness(X_dry, X_est);        
            scores.CC(rv,jk,n)    = gof.CC;
            scores.mse(rv,jk,n)   = gof.mse;
            scores.nmse(rv,jk,n)  = gof.nmse;
            scores.CCf(rv,jk,n,:) = gof.CCf;            
            sti.est(rv,jk,n)      = STI(X_est, spec_st);
            sti.est_probe(rv,jk,n)= STI({X_dry, X_est}, spec_st);
            
            % For comparing between stimuli of various DRR conditions and the
            % reconstructed spectrograms
            X_drr                 = spec_st.Sft{rv}(:, jk == splits.idx);            
            gof2                  = goodness(X_drr, X_est);        
            scores2.CC(rv,jk,n)   = gof2.CC;
            scores2.mse(rv,jk,n)  = gof2.mse;
            scores2.nmse(rv,jk,n) = gof2.nmse;
            scores2.CCf(rv,jk,n,:)= gof2.CCf;
            %sti2.est(rv,jk,n)     = STI({X_est, X_drr}, spec_st);
            
            % * # of neurons doesn't make any difference for the STIMULUS
            if 1 == n
                sti.test(rv,jk) = STI(X_drr, spec_st);
                sti.test_probe(rv,jk) = STI({X_drr, X_dry}, spec_st);
            end
            
        end
    end
    

end
 


%% Statistics
% Dims: [drr x # splits x # units]
scores.mu.CC = squeeze( median(scores.CC, 2) );
scores.SE.CC = squeeze( mad(scores.CC,[],2)/sqrt(n_JK) );
scores.mu.mse = squeeze( median(scores.mse, 2) );
scores.SE.mse = squeeze( mad(scores.mse,[],2)/sqrt(n_JK) );
scores.mu.nmse = squeeze( median(scores.nmse, 2) );
scores.SE.nmse = squeeze( mad(scores.nmse,[],2)/sqrt(n_JK) );

scores2.mu.CC = squeeze( median(scores2.CC, 2) );
scores2.SE.CC = squeeze( mad(scores2.CC,[],2)/sqrt(n_JK) );
scores2.mu.mse = squeeze( median(scores2.mse, 2) );
scores2.SE.mse = squeeze( mad(scores2.mse,[],2)/sqrt(n_JK) );
scores2.mu.nmse = squeeze( median(scores2.nmse, 2) );
scores2.SE.nmse = squeeze( mad(scores2.nmse,[],2)/sqrt(n_JK) );



%% Save for later use
%{
analyze.path = '../_data/Analyze/';
analyze.filename = sprintf('Analysis_%s_(%s)_units(%d)_bw(%g)_fbands(%d)',...
    finfo.type, date, n_neurons, binwidth, n_bands);
save([analyze.path, analyze.filename], '-v7.3',...
    'binwidth', 'scores', 'scores2',...
    ... 'obj_list',...
    'spec_st', 'H', 'splits', 'n_neurons', 'CC_broadband_envelope',...
    'tbl_impale', 'tbl_metadata');
%}


%% Save data for the STI
%{
analyze.path = '../_data/STI/';
analyze.filename = sprintf('STI_%s_(%s)_units(%d)_bw(%g)_fbands(%d)',...
    finfo.type, date, n_neurons, binwidth, n_bands);
save([analyze.path, analyze.filename],...
    'sti', ...
    'binwidth', 'spec_st', 'H', 'splits',...
    'tbl_impale', 'tbl_metadata');
%}


%% Correlation Coefficients between STIMULI
%{
Sdry = spec_st.Sft{drr.dry};
CCs = nan(1, n_drr);    % CCs of responses
CCs_std = nan(1, n_drr);    % CCs of responses
for k = 1:n_drr
    Sk = spec_st.Sft{k};
    
    % CCs: correlation between SRY & DRR stimuli
    cs_k = diag( zca(Sdry')' * zca(Sk') );
    CCs(k) = mean(cs_k)/size(Sk,2);
    CCs_std(k) = sqrt(mean((cs_k/size(Sk,2) - CCs(k)).^2));
    
end
%}

Sdry = spec_st.Sft{drr.dry};
CCs = nan(1, n_drr);    % CCs of responses
for k = 1:n_drr
    Sk = spec_st.Sft{k};
    
    % CCs: correlation between SRY & DRR stimuli    
    gof = goodness(Sdry, Sk);
    CCs(k) = gof.CC;
end




%% Statistics
figure(0+fignum);
clf;

fontsize = 12;
markersize = 34;

x = 1:n_drr;
M = [scores.mu.CC(drr.ordered,end), scores2.mu.CC(drr.ordered,end)];
errbar = [scores.SE.CC(drr.ordered,end), scores2.SE.CC(drr.ordered,end)];
h = bar(M); 
set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.ordered));
h(1).FaceColor = aux.rpalette(sprintf('new%02d',1));
h(2).FaceColor = aux.rpalette(sprintf('new%02d',2));

dx = 0.2*h(1).BarWidth;
hold on
% ERROR BARs
er = errorbar([(1:5)'-dx, (1:5)'+dx], M, errbar);

% CCs
plth = plot(1:n_drr, CCs(drr.ordered), 'ks:',...
    'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k' );

hold off
for k = 1:length(er)
    er(k).Color = [0 0 0];                            
    er(k).LineStyle = 'none'; 
    er(k).LineWidth = 2;
end
ylim([0.0, 1.0]);
legend(h, {'To dry', 'To DRR'}, 'Location', 'northeastoutside');
ylabel('CC');
xlabel('DRR');
aux.ctitle( sprintf('Comparing Reconstructed Spectrograms of %s', finfo.type),...
    'to Dry and Other DRR Conditioned Stimuli');
set(gca, 'FontSize', fontsize);



%% Correlation Coefficients between STIMULI
figure(1+fignum);
clf;

fontsize = 12;
markersize = 34;

% Sdry = spec_st.Sft{drr.dry};
% CCs = nan(1, n_drr);    % CCs of responses
% CCs_std = nan(1, n_drr);    % CCs of responses
% for k = 1:n_drr
%     Sk = spec_st.Sft{k};
%     
%     % CCs: correlation between SRY & DRR stimuli
%     cs_k = diag( zca(Sdry')' * zca(Sk') );
%     CCs(k) = mean(cs_k)/size(Sk,2);
%     CCs_std(k) = sqrt(mean((cs_k/size(Sk,2) - CCs(k)).^2));
%     
% end

% plth = plot(1:n_drr, CCs(drr.ordered), 'ks:',...
%     'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k' );
% hold on
% CC of broadband envelopes (DRY to DRR)
% plth(2) = plot(CC_broadband_envelope, 'dk:', 'MarkerSize', 0.4*markersize);

% ERRORBAR
CC_stim_ordered = CCs(drr.ordered)';
dCCs = CC_stim_ordered./CC_stim_ordered(1);
dM =  M./M(1,:);

plth = plot(1:drr.n_drr, dCCs, 'ks:', 'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k');
hold on
plot(repmat((1:drr.n_drr)', 1,2) + 0.05*[-1, 1], dM, '.');
hold off


% ylim([0.5, 1.2]);
set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.ordered));
ylabel('CC');
xlabel('DRR');
% legend([ehnd, plth],...
%     {'To Dry (norm. by dry)', 'To DRR (norm. by dry)', ...
%     '$\overline{CC}_{(dry-to-DRR)}$ (spectrogram envelopes)', ...
%     '$CC_{(dry-to-DRR)}$ (broadband envelope)'}, ...
%     'Location', 'southwest');

set(gca, 'FontSize', fontsize);



%% Plot: comparing between DRY & other DRR reconstructions 
I = 2;                  % split number\jackknife number to use
rvi = drr.sortby(1);    % reverberant condition
X_test = spec_st.Sft{drr.dry}(:, I == splits.idx);          % DRY

% Whitening 
X_test_ = zca( X_test' )';

fkHz= 1e-3*spec_st.f;                    % (kHz)
t   = spec_st.t(1:size(X_test,2));       % (sec)
f_slc = [10, 15, 28];   % index number


% Selected frequency band to show
fprintf('\n===================\n');
fprintf('--> n_neurons  : %d\n',  obj_list{rvi,I}.n_neurons);


% fprintf('\n-> Envelope #%d\n', q);
cc_env = nan(n_drr, 3);
for k = 1:n_drr
    rvi = drr.sortby(k);    % reverberant condition
    X_drr  = spec_st.Sft{rvi}(:, I == splits.idx);    % -8.2 dB; extreme reverberation
    X_est  = obj_list{rvi,I}.X_est;
    X_drr_  = zca( X_drr' )';
    X_est_  = zca( X_est' )';

    cc_env(k,1) = mean(diag(X_test_*X_est_'))/size(X_test_,2);
    cc_env(k,2) = mean(diag(X_test_*X_drr_'))/size(X_test_,2);
    cc_env(k,3) = mean(diag(X_drr_*X_est_'))/size(X_test_,2);

end
fprintf('--> cc_env\n');
tbl_cc_env = array2table(cc_env,...
    'VariableNames', {'dry2est', 'dry2drr', 'drr2est'},...
    'RowNames', drr.labels(drr.sortby(1:n_drr)));
disp(tbl_cc_env);
% writetable(tbl_cc_env, '../_data/tbl_cc_env.xls', 'WriteVariableNames', 1, 'WriteRowNames', 1);



%% CCs & NMSEs
figure(3+fignum);
clf;
fontsize = 24;
markersize = 20;

drr_idx = drr.sortby(1:n_drr);

% Automatically choose the data with most (last) of the neurons
n_neurons_idx = (finfo.n_neurons == n_neurons); 
n_units = size(H,3);
assert( scores.n_units(n_neurons_idx) == n_units, '--> WRONG # of units!!!');


% CC
ax = subplot(1,2,1);
D = scores.CC(drr_idx,:,n_neurons_idx)';
plth = plot(D, 's:',...
    'MarkerSize', markersize);
if strcmpi('MUA', finfo.type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end
xlabel('Speaker Sex');
ylabel('CC');
% ylabel('$CC(S_{DRY}(t),\hat{S}_{DRR}(t))$');
% title(sprintf('$CC_{RS}^{%s}$ (%d %ss)', finfo.type, n_units, finfo.type));
% title(sprintf('$CC(S_{DRY}(t),\\hat{S}_{DRR}(t))$ (%d %ss)', n_units, finfo.type));
aux.ctitle('Dry Stimulus vs. Reconstruction',...
    sprintf('(%d %ss)', n_units, finfo.type));

set(gca, 'FontSize', fontsize);
legend( drr.labels{drr_idx}, 'Location', 'southeast');
axis tight 
ylim([0, 1]);

speaker_sex_list = {'M$_1$','M$_2$','M$_3$','F$_1$','M$_4$','F$_2$','F$_3$',...
    'F$_4$','M$_5$','F$_5$','F$_6$','M$_6$'};
xticks = get(ax(1), 'XTick');
set(ax(1), 'XTickLabel', speaker_sex_list(xticks));

% NMSE
ax(2) = subplot(1,2,2);
plth = plot(scores.nmse(drr_idx,:,n_neurons_idx)', 's:', 'MarkerSize', markersize);
if strcmpi('MUA', finfo.type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end
xlabel('Speaker Sex');
ylabel('NMSE');
title(sprintf('$NMSE_{RS}^{%s}$', finfo.type));
set(gca, 'FontSize', fontsize);
axis tight 
ylim([0, 1]);
xticks = get(ax(2), 'XTick');
set(ax(2), 'XTickLabel', speaker_sex_list(xticks));

aux.abc(ax);





%% Plot: Compare CCs vs. Speaker number
% Correlations between responses to different DRRs
figure(5+fignum);
clf;
fontsize = 30;
markersize = 40;

% n_drr = 5;
drr_idx = drr.sortby(1:n_drr);
% Automatically choose the data with most (last) of the neurons
n_neurons_idx = (finfo.n_neurons == n_neurons);     

CCs = nan(1, n_drr);    % CCs of stimuli
CCr = nan(1, n_drr);    % CCs of responses

% Use only "valid" measurements, i.e., measurements that include all 5 DRR
% conditions
valid_idx = logical( squeeze(prod(~isnan(sum(H(:,1:n_drr,:),1)),2)) );
n_units = size(H,3);
assert(n_units == nnz(valid_idx), '--> Some of the loaded units are INVALIDE...');


Sdry = spec_st.Sft{drr.dry};
Hdry = squeeze(zca(H(:,drr.dry,valid_idx)))';
for k = 1:n_drr
    Sk = spec_st.Sft{k};
    
    % CCs: correlation between SRY & DRR stimuli
    cs_k = diag( zca(Sdry')' * zca(Sk') );
    CCs(k) = mean(cs_k)/size(Sk,2);
    
    cr_k = diag( Hdry * squeeze(zca(H(:,k,valid_idx))) );
    CCr(k) = mean(cr_k)/(size(H,1));
    
end

% CC(DRY, DRR)
plth = plot(1:n_drr, CCs(drr_idx), 'sk:', 'MarkerSize', 0.4*markersize);
arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));

hold on
% CC of broadband envelopes (DRY to DRR)
plot(CC_broadband_envelope, 'dk:', 'MarkerSize', 0.4*markersize);

% Response CCr(dry, DRR)
plth = plot(1:n_drr, CCr(drr_idx), 's:', 'MarkerSize', 0.4*markersize);
if strcmp('MUA', finfo.type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end

% Mean CC reconstructions
plth = plot(scores.mu.CC(drr_idx, n_neurons_idx), 's:', 'MarkerSize', 0.4*markersize);
if strcmp('MUA', finfo.type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end

hold off

set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.sortby));
xlabel('DRR');
ylabel('CC');
title(sprintf('Compare CC (%d %ss)', n_neurons, finfo.type));
ylim([0.0, 1.0]);
legend('$CC_{stimulus}(S_{dry}(f,t), S_{DRR}(f,t))$', ...
    'CC (Stim BB Envelope)', ...    
    sprintf('$CC_{response}(r_{dry}, r_{DRR})$ (%d %ss)',n_neurons,  finfo.type),...
    sprintf('$CC_{RS}(dry, DRR)$ (%d %ss)',n_neurons,  finfo.type),...
    'Location', 'southwest');
set(gca, 'FontSize', fontsize);






%% Plot: CC & NMSE vs. DRRs
% CC:
figure(7+fignum);
clf;

fontsize = 28;
fs_legend= 20;
markersize = 18;

drr_idx = drr.sortby(1:n_drr);

% CC
ax = subplot(1,2,1);
plth = plot(scores.mu.CC(drr_idx,:), 's:', 'MarkerSize', markersize);
if strcmpi('MUA', finfo.type)   
    colors = get(plth,'Color');
    if 1 < size(colors,1)
        arrayfun(@(I) set(plth(I), 'MarkerFaceColor', colors{I}), 1:length(colors) );
    else
        set(plth, 'MarkerFaceColor', colors);
    end
end
hold on
plot(CC_broadband_envelope, 'dk:', 'MarkerSize', markersize);
plot(CCs(drr_idx), 'sk:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k');
hold off

set(gca, 'XTickLabel', drr.labels(drr.sortby));
xlabel('DRR');
ylabel('CC($S_{dry}(f,t)$-to-$\hat{S}(f,t)$)');
title(sprintf('$CC_{RS}^{%s}$', finfo.type));
set(gca, 'FontSize', fontsize);
units_list = arrayfun(@(n) sprintf('%d %ss', n,  finfo.type), finfo.n_neurons, 'UniformOutput', 0);
legend_list = [units_list, 'CC (Stim BB Envelope)', 'CC (Stim Dry-to-DRR)'];
legend_h = legend(legend_list, 'Location', 'southwest');
legend_h.FontSize = fs_legend;

ylim([0.4, 1.0]);

% NMSE
ax(2) = subplot(1,2,2);
plot(scores.mu.nmse(drr_idx,:), '.:', 'MarkerSize', markersize);
set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.sortby));
xlabel('DRR');
ylabel('$NMSE_{(S_{(f,t)}, \hat{S}_{(f,t)})}$');
title(sprintf('$NMSE_{RS}^{%s}$', finfo.type));
set(gca, 'FontSize', fontsize);
aux.abc(ax);
axis tight



%% BoxPlot: CC vs. DRRs
% CC:
figure(9+fignum);
clf;

fontsize = 28;
fs_legend= 20;
markersize = 18;

drr_idx = drr.sortby(1:n_drr);

% CC
ax = subplot(1,2,2);
% plth = plot(scores.mu.CC(drr_idx,:), 's:', 'MarkerSize', markersize);
plth = boxplot(scores.CC(drr_idx,:)', 'Notch', 'off', 'colors', aux.rpalette(4));
if strcmpi('MUA', finfo.type)   
    colors = get(plth,'Color');
    if 1 < size(colors,1)
        arrayfun(@(I) set(plth(I), 'MarkerFaceColor', colors{I}), 1:length(colors) );
    else
        set(plth, 'MarkerFaceColor', colors);
    end
end
hold on
% plot(CC_broadband_envelope, 'dk:', 'MarkerSize', markersize);
plot(CCs(drr_idx), 'sk:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k');
hold off

set(gca, 'XTickLabel', drr.labels(drr.sortby));
xlabel('DRR');
ylabel('CC($\hat{S}_{dry}(f,t)$-to-$\hat{S}(f,t)$)');
% title(sprintf('$CC_{RS}^{%s}$', finfo.type));
set(gca, 'FontSize', fontsize);
units_list = arrayfun(@(n) sprintf('%d %ss', n,  finfo.type), finfo.n_neurons, 'UniformOutput', 0);
legend_list = [units_list, 'CC (Stim BB Envelope)', 'CC (Stim Dry-to-DRR)'];
legend_h = legend(legend_list, 'Location', 'southwest');
legend_h.FontSize = fs_legend;

ylim( [0.4, 1.0] );

% % NMSE
% ax(2) = subplot(1,2,2);
% plot(scores.mu.nmse(drr_idx,:), '.:', 'MarkerSize', markersize);
% set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.sortby));
% xlabel('DRR');
% ylabel('$NMSE_{(S_{(f,t)}, \hat{S}_{(f,t)})}$');
% title(sprintf('$NMSE_{RS}^{%s}$', finfo.type));
% set(gca, 'FontSize', fontsize);
% aux.abc(ax);
% axis tight




%% Correlations between frequencies of selected intervals (jackknife)
figure(10+fignum);
clf;

fontsize = 32;
markersize = 30;

% chunk_number  = 12;
drr_idx       = drr.sortby(1:n_drr);
n_neurons_idx = find(finfo.n_neurons == n_neurons);

CCf_intervals = squeeze( scores.CCf(drr_idx(1),:,n_neurons_idx,:) );

% CC Averaged Across Intervals
% ax(2) = subplot(1,3,1);
ax = gca;
slc_interval = [2, 5, 9, 12];   % selected inteval number to show
plth = plot(1e-3*f, CCf_intervals(slc_interval,:)', '.:', 'MarkerSize', markersize);
xlabel('Frequency (kHz)');
ylabel('$CC_{(S_{dry}(f,t), \hat{S}_{dry}(f,t))}$');
title(sprintf('$CC_{RS}^{%s}$ Across Intervals (%d %ss)', finfo.type, n_neurons, finfo.type));
set(gca, 'FontSize', fontsize);
set(ax(1), 'FontSize', fontsize);

hold on
mu    = median(CCf_intervals,1);
sigma = mad(CCf_intervals);
plth(end+1) = plot(1e-3*f, mu, ':k');
yarea = [mu(:) - sigma(:), 2*sigma(:)];
harea = area(1e-3*f, yarea);
hold off
harea = harea(end:-1:1);
harea(2).Visible = 'off';
harea(1).FaceAlpha = 0.4;
harea(1).EdgeAlpha = 0.4;
harea(1).FaceColor = 'k';

legend_str = arrayfun(@(X) sprintf('Jackknife Speaker %d', X), slc_interval, 'UniformOutput', 0);
legend_str = [legend_str, 'median', 'mean absolute deviation'];
legend([plth; harea(1)], legend_str, 'Location', 'southeast', 'FontSize', 0.75*fontsize);

ylim( [0.7, 1.0] );
xlim( [0.85*1e-3*f(1), 1.01*1e-3*f(end)] );


%% CC (DRY-to-DRR) vs. CC (DRR-to-DRR)
figure(12+fignum);
clf;

fontsize = 24;
markersize = 14;

CCf_avg_intervals = squeeze( mean(scores.CCf(drr_idx,:,n_neurons_idx,:), 2) );
CCf2_avg_intervals = squeeze( mean(scores2.CCf(drr_idx,:,n_neurons_idx,:), 2) );

% CCf (DRY vs. estimated-DRR) Averaged Across Intervals
ax = subplot(1,2,2);
plth = plot(1e-3*f, CCf_avg_intervals', 's:', 'MarkerSize', markersize);
if strcmpi('MUA', finfo.type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end

% axis tight
xlim([0, 7.6]);
xlabel('Frequency (kHz)');
title( sprintf('$CC^{%s}_{RS}$ Dry-vs-DRR', finfo.type) );
ylabel('$CC_{(S_{\bf{dry}}(f,t), \hat{S}_{DRR}(f,t))}$');

set(gca, 'FontSize', fontsize);

% CCf2 (DRR vs. estimated-DRR) Averaged Across Intervals
ax(2) = subplot(1,2,1);
plth = plot(1e-3*f, CCf2_avg_intervals', 's:', 'MarkerSize', markersize);
if strcmpi('MUA', finfo.type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end

xlim([0, 7.6]);
xlabel('Frequency (kHz)');
ylabel('$CC_{(S_{\bf{DRR}}(f,t), \hat{S}_{DRR}(f,t))}$');
title(sprintf('$CC_{RS}^{%s}$ DRR-vs-DRR', finfo.type));

set(gca, 'FontSize', fontsize);
legend(plth, drr.labels(drr_idx), 'Location', 'southeast', 'FontSize', fontsize);
aux.abc(ax(end:-1:1));

linkaxes(ax);




%% Plot: spectrograms of stimuli vs. reconstructions
figh = figure(20+fignum);
clf;

spec_thr    = 20;  % (dB) threshold\remove from spectrogram before plotting it
color_axis  = [0, 45]; 
fontsize    = 24;

I = 6;      % split number \ jackknife number to use
nt      = nnz(I == splits.idx);
t   	= spec_st.t(1:nt);       % (sec)
X_test_ = nan(n_bands, nt);
X_est_  = nan(n_bands, nt);

for rv = 1:n_drr
    X_test = spec_st.Sft{rv}(:, I == splits.idx);   
    X_test_(:,:,rv) = max(0, X_test - spec_thr);    
    %X_test_(:,:,rv) = zca(X_test')';
    
    % Whitening 
    X_est  = obj_list{rv,I}.X_est;
    X_est_(:,:,rv) = max(0, X_est - spec_thr);
    %X_est_(:,:,rv) = zca(X_est')';    
end
            
ax = zeros(5,2);
% fs = 1/(1e-3*binwidth);
% t  = spec_st.t(1:size(X_test,2)); %1/fs*[1:size(X_test,2)]';
% f  = 1e-3*spec_st.f;

for rv = 1:n_drr    
    rvi = drr.sortby(rv);

    %ax(rv,1) = subplot(5,2,1+2*(rv-1));            
    ax(rv,1) = subplot(5,2, 5*2 - (1+2*(rv-1)));    
    spec.plot_spectrogram(t, 1e-3*f, X_test_(:,:,rvi), figh, 1);
    
    %ax(rv,2) = subplot(5,2,2*rv);
    ax(rv,2) = subplot(5,2, 5*2-2*(rv-1));
    spec.plot_spectrogram(t, 1e-3*f, X_est_(:,:,rvi), figh, 1);
    
    % Add the DRRs to the ylabels
    if 3 ~= rv
        ylabel(ax(rv,1), drr.labels{rvi});
    else
        ystr = aux.ctitle('Frequency (kHz)\\', drr.labels{rvi});
        ylabel(ax(rv,1), ystr);        
    end

end

ax = flipud(ax);
arrayfun(@(X) caxis(X, color_axis), ax);   % (dB) set the range of colors

set(ax, 'FontSize', fontsize);
set(ax(1:4,:), 'XTickLabel', '');
set(ax(:,2), 'YTickLabel', '');
xlabel(ax(end,1), 'Time (sec)');
xlabel(ax(end,2), 'Time (sec)');
title(ax(1,1), 'Stimulus');
title(ax(1,2), sprintf('Reconstruction (%d %ss, dry filters)', n_units, finfo.type));
aux.abc(ax(1,:), 'fontsize', 2*fontsize, 'location', 'northwestoutside');
linkaxes(ax);       
            
% Add a colorbar and make sure that the axes doesn't change in size (push 
% the colorbar "outside").
if 1
    ax_j = ax(end,end);
    pos = get(ax_j, 'Position');
    
    colorbar(ax_j);
    set(ax_j, 'Position', pos);
    
end
     



    

%% Plot: RMD, MG, Kurtosis, and CCr
% RMD: response modulation depth
% MG : modulation gain
% CCr: correlation coefficient between DRY and other DRR conditions
%n_valid = nnz(valid_idx);   % # of valid measurements
fontsize = 32;
markersize = 30;

n_units = size(H,3);

drr_idx = drr.sortby(1:n_drr);
% n_neurons_idx = find(finfo.n_neurons == n_neurons);
n_drr = 5;

% Select data for analysis
D           = squeeze( H(:,:, 1:n_units) );
scores.CCr  = nan(n_units, n_drr);       % % CC(response dry & response)
scores.MDr  = nan(n_units, n_drr);
scores.MG   = nan(n_units, n_drr);
scores.skew = nan(n_units, n_drr);
scores.ku   = nan(n_units, n_drr);

% if strcmpi('SU', finfo.type)
%     % Select only measurements with single units (SUs)
%     tbl_36sec = tbl_impale(1 == tbl_impale.SPK & tbl_impale.duration_sec == duration_sec, :);
%     units_all_valid = find(data.slc.valid_neurons == 1);
%     units_used = units_all_valid(1:n_units);
%     tbl_valid = tbl_36sec(units_used, :);
% else
%     tbl_36sec = tbl_impale(tbl_impale.duration_sec == duration_sec, :);
%     units_all_valid = find(data.slc.valid_neurons == 1);
%     units_used = units_all_valid(1:n_units);
%     tbl_valid = tbl_36sec(units_used, :);
% end


tbl_valid = tbl_impale( data.slc.unit_used, : );
assert( 0==nnz(tbl_valid.CF==0) && 0 == nnz(isnan(tbl_valid.CF)) );

for n = 1:n_units
    % Option #1:
    % Stimulus envelope; find the closest frequency band to the neuron's CF 
    [~, idx_cf] = min(abs(f - tbl_valid.CF(n)));

     for rv = 1:n_drr
         % CC
        dummy = corrcoef(D(:,drr.dry,n), D(:,rv,n));
        scores.CCr(n,rv) = dummy(1,2);
        
        % Response MD (modulation depth)
        scores.RMD(n,rv) = MD(D(:,rv,n), 1);
        
        % Envelope's MD 
        % Option #1: using the spectrogram
        % %{
        env_db = spec_st.Sft{rv}(idx_cf,:);        
        cfloor = spec_st.db_floor;              % The spectrogram (Sft) amplitudes 
        cmax = spec_st.max_Sft;                 % are of log scale, so we need to 
        env_ = cmax*10.^((env_db + cfloor)/20); % revert it
        MD_env = MD( env_ );
        %}
        % Option #2: using the stimulus (more computation time)
        %{
        stim_env_rv = calc_stimulus_envelope(stim_st.Y(:,rv), tbl_valid.CF(n), fs, stim_st.fs);
        MD_env = MD( stim_env_rv, 1 );
        %}
        
        % MG (modulation gain) = RMD/MDenv
        scores.MG(n,rv) = db(scores.RMD(n,rv)./MD_env);
        scores.skew(n,rv) = skewness( D(:,rv,n) );
        scores.kr(n,rv) = kurtosis( D(:,rv,n) );

%                 '### LPC DEBUG ###' 
%                 x = env_;
%                 a = lpc(x, 3);
%                 est_x = filter([0 -a(2:end)], 1, x);
%                 e = x-est_x;
%                 corrcoef(e(:), D(:,3,n)-mean(D(:,3,n)))
%                 plot([e(:), (D(:,3,n)-mean(D(:,3,n)))])
     end
end

% % Fisher transformation
% scores.RMD = tanh(scores.RMD);
% scores.MG = tanh(scores.MG);
% scores.kr = tanh(scores.kr);
% scores.CCr = tanh(scores.CCr);


%
figure(30+fignum);
clf;

fontsize = 20;

ax = subplot(4,1,1);
plot_dotbox(scores.RMD(:,drr_idx), 'labels', drr.labels(drr_idx));
ylabel(aux.ctitle('RMD', '$(\sqrt{2}\sigma/\mu)$'));
set(gca, 'XTickLabel', '');
xlabel('');
title(ax(1), sprintf('%s', finfo.type));

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
plot_dotbox(scores.CCr(:,drr_idx), 'labels', drr.labels(drr_idx));
ylabel('$CC$');

hold on
plth = plot(CC_broadband_envelope, 'dk:', 'MarkerSize', 0.4*markersize);
plth(2) = plot(CCs(drr_idx), 'sk:', 'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k');
hold off
legend(plth, {'CCe broadband envelope', 'CCs stimulus'});

xlabel('Direct to Reverberation Ratio (dB)');
set(ax, 'FontSize', fontsize);

aux.abc(ax, 'fontsize', 2*fontsize, 'location', 'northwestoutside');


%% Friedman statistical test
clear pvalues
% alpha = 0.05;
display = 'off';
fontsize = 32;

% %{
figure(35+fignum);
clf;
ax = subplot(2,2,1);
[p,tbl,stats] = friedman(scores.RMD(:,drr_idx),1,'off');
[c, pvalues.RMD] = multcompare(stats, 'Display', display);
errorbar(pvalues.RMD(:,1), 1:n_drr, pvalues.RMD(:,2),...
    'o', 'horizontal', 'LineWidth', 2);
title(sprintf('RMD (%s)', finfo.type));
set(gca, 'YTickLabel', arrayfun(@(I) drr.labels{I}, drr_idx, 'uni', 0)' );

ax(2) = subplot(2,2,2);
[p,tbl,stats] = friedman(scores.MG(:,drr_idx),1,'off');
[c, pvalues.MG, ~, ~] = multcompare(stats, 'Display', display);
errorbar(pvalues.MG(:,1), 1:n_drr, pvalues.MG(:,2),...
    'o', 'horizontal', 'LineWidth', 2);
title(sprintf('MG (%s)', finfo.type));
% set(gca, 'YTickLabel', '');
set(gca, 'YTickLabel', arrayfun(@(I) drr.labels{I}, drr_idx, 'uni', 0)' );

ax(3) = subplot(2,2,3);
[p,tbl,stats] = friedman(scores.kr(:,drr_idx),1,'off');
[c, pvalues.kr, ~, ~] = multcompare(stats, 'Display', display);
errorbar(pvalues.kr(:,1), 1:n_drr, pvalues.kr(:,2),...
    'o', 'horizontal', 'LineWidth', 2);
title(sprintf('Kurtosis (%s)', finfo.type));
set(gca, 'YTickLabel', arrayfun(@(I) drr.labels{I}, drr_idx, 'uni', 0)' );

ax(4) = subplot(2,2,4);
[p,tbl,stats] = friedman(scores.CCr(:,drr_idx),1,'off');
[c, pvalues.CC, ~, ~] = multcompare(stats, 'Display', display);
errorbar(pvalues.CC(:,1), 1:n_drr, pvalues.CC(:,2),...
    'o', 'horizontal', 'LineWidth', 2);
title(sprintf('CC (%s)', finfo.type));
% set(gca, 'YTickLabel', '');
set(gca, 'YTickLabel', arrayfun(@(I) drr.labels{I}, drr_idx, 'uni', 0)' );
xlabel('Mean column rank');

set(ax, 'Fontsize', 0.8*fontsize);
axis tight
aux.abc(ax, 'fontsize', 1.5*fontsize);
%}

%{
pvalues.RMD= multiple_friedman_test(scores.RMD(:,drr_idx));
pvalues.MG = multiple_friedman_test(scores.MG(:,drr_idx));
pvalues.kr = multiple_friedman_test(scores.kr(:,drr_idx));
pvalues.CC = multiple_friedman_test(scores.CC(:,drr_idx));

pvalues.tbl_all = [pvalues.RMD; pvalues.MG; pvalues.kr; pvalues.CC];
pvalues.tbl_all.Properties.RowNames = {'RMD', 'MG', 'Kr', 'CC'};

fprintf('-> Friedman statistical test (%s)\n', finfo.type);
disp(pvalues.tbl_all);


writetable(pvalues.tbl_all, sprintf('multiple_friedman_test_%s.xls', finfo.type),'WriteRowNames',true);
%}


%% Plot: STI
figure(40+fignum);
clf;

fontsize = 24;
markersize = 24;

drr_idx = drr.sortby(1:n_drr);
n_neurons_idx = find(finfo.n_neurons == n_neurons);
n_drr = 5;

ax = subplot(2,1,1);
plth = plot(sti.test(drr_idx,:)', 'x', 'MarkerSize', 0.75*markersize);
xlabel('Speaker Number');
ylabel('STI');
hold on
plth2 = plot(sti.est(drr_idx,:,n_neurons_idx)', '.:', 'MarkerSize', markersize);
hold off
for q = 1:n_drr
    plth2(q).Color = plth(q).Color;
end
legend([plth2; plth(1)], [drr.labels(drr_idx), 'Test Stimulus'], 'Location', 'northeastoutside');
title(sprintf('Speech Transmission Index (%d %ss)', n_neurons, finfo.type));
xlim([0, 13]);

ax(2) = subplot(2,1,2);
plth2 = plot(squeeze(mean(sti.est(drr_idx,:,:),2)), 's:', 'MarkerSize', 0.5*markersize);
plth2(end).MarkerFaceColor = plth2(end).Color;  
set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.sortby));
xlabel('DRR');
xlim([0.75, 1.05*n_drr]);
ylabel('STI');
hold on
plth = plot(mean(sti.test(drr_idx,:),2), 'x',...
    'MarkerSize', markersize, 'LineWidth', 3, 'Color', aux.rpalette('new01'));
hold off
units_str = arrayfun(@(I) sprintf('%d meas.',finfo.n_neurons(I)), 1:n_neuron_list,'UniformOutput',0);
legend([plth2; plth], [units_str, 'Test Stimulus'], 'Location', 'northeastoutside');

set(ax, 'FontSize', fontsize);
aux.abc(ax);





%% Plot: STI with a probe
figure(45+fignum);
clf;

fontsize = 24;
markersize = 24;

drr_idx = drr.sortby(1:n_drr);
n_neurons_idx = find(finfo.n_neurons == n_neurons);
n_drr = 5;

ax = subplot(2,1,1);
plth = plot(sti.test_probe(drr_idx,:)', 'x', 'MarkerSize', 0.75*markersize);
xlabel('Speaker Number');
ylabel('STI');
hold on
plth2 = plot(sti.est_probe(drr_idx,:,n_neurons_idx)', '.:', 'MarkerSize', markersize);
hold off
for q = 1:n_drr
    plth2(q).Color = plth(q).Color;
end
legend([plth2; plth(1)], [drr.labels(drr_idx), 'Test Stimulus'], 'Location', 'northeastoutside');
title(sprintf('Speech Transmission Index (%d %ss)', n_neurons, finfo.type));
xlim([0, 13]);

ax(2) = subplot(2,1,2);
plth2 = plot(squeeze(mean(sti.est_probe(drr_idx,:,:),2)), 's:', 'MarkerSize', 0.5*markersize);
plth2(end).MarkerFaceColor = plth2(end).Color;  
set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.sortby));
xlabel('DRR');
xlim([0.75, 1.05*n_drr]);
ylabel('STI');
hold on
plth = plot(mean(sti.test_probe(drr_idx,:),2), 'x',...
    'MarkerSize', markersize, 'LineWidth', 3, 'Color', aux.rpalette('new01'));
hold off
units_str = arrayfun(@(I) sprintf('%d meas.',finfo.n_neurons(I)), 1:n_neuron_list,'UniformOutput',0);
legend([plth2; plth], [units_str, 'Test Stimulus'], 'Location', 'northeastoutside');

set(ax, 'FontSize', fontsize);
aux.abc(ax);





%% Plot: Speech-envelope spectrum & Modulation index
figh = figure(52+fignum);
clf;

fontsize = 32;

chunk_number = 5;
% chunk_number= 1:12,     
% chunk_number= [1, 2, 3, 5, 9, 12];                      % Male speakers
% chunk_number= setdiff(1:n_splits, chunk_number); % Female speakers
fun_stim    = @(n,k) spec_st.Sft{n}(:, k == splits.idx);
fun_drr     = @(n,k) obj_list{n,k}.X_est;
n_chunks    = length(chunk_number);
n_drr       = 5;
sub_y       = 2;
sub_x       = n_drr;

if ~exist('tbl_metadata', 'var')
    dummy = load('C:\Users\barzelo\Google Drive\codeOnCloud\Reverberation\Stimulus\.database\Project1\metadata_(36)_wav_(30-Jun-2020).mat');
    tbl_metadata = dummy.tbl_metadata;
end
%   tbl_metadata(:,'txt')
disp([tbl_metadata.txt{chunk_number}]);


ax = nan(sub_y, sub_x);
for rv = 1:n_drr      
    rvi = drr_idx(rv);
    
    Xs = [];
    Xdrr = [];

    % Aggregate chunks
    for cn = 1:n_chunks
        cn_idx = chunk_number(cn);
        
        win = 1;
        %win = hanning(size(Xs_i,2))';
        %win = tukeywin(size(Xs_i,2), 0.1)';
        
        Xs_i = fun_stim(rvi,cn_idx);        
        Xs   = [Xs  , Xs_i.*win];
        
        Xdrr_i = fun_drr(rvi,cn_idx);        
        Xdrr = [Xdrr, Xdrr_i.*win];
    end
    
    % MDT & STI        
    % Stimulus modulation index
    [~, MIs, fm] = STI(Xs, spec_st);
    %[~, MIs, fm] = STI(Xs, binwidth);

    % Reconstructed stimulus (RS) modulation index
    [~, MIdrr, ~] = STI(Xdrr, spec_st);
    %[~, MIdrr, ~] = STI(Xdrr, binwidth);

    I = 1;
    sub_idx = (I-1)*n_drr + rv;
    ax(I,rv) = subplot(sub_y, sub_x, sub_idx);
    spec.plot_spectrogram(fm, 1e-3*f, MIs, figh, [], 0.75*fontsize);
    colorbar(ax(I, rv), 'off');
    xlim([fm(1), fm(66)]);  % (Hz)

    sub2_idx = I*n_drr + rv;
    ax(1+I,rv) = subplot(sub_y, sub_x, sub2_idx);
    spec.plot_spectrogram(fm, 1e-3*f, MIdrr, figh, [], 0.75*fontsize);
    colorbar(ax(1+I, rv), 'off');
    xlim([fm(1), fm(66)]);  % (Hz)

    if 1 < rv 
        set(ax(I,rv),'XTickLabel', '');            
        set(ax(1+I,rv),'XTickLabel', '');   
        xlabel(ax(I,rv), '');
        xlabel(ax(1+I,rv), '');

        set(ax(I,rv),'YTickLabel', '');
        set(ax(1+I,rv),'YTickLabel', '');
        ylabel(ax(I,rv), '');
        ylabel(ax(1+I,rv), '');
    else
        set(ax(I,rv),'XTickLabel', '');            
        xlabel(ax(I,rv), '');
        set(ax(I,rv),'YTickLabel', '');
        ylabel(ax(I,rv), '');            
    end

    if 1 == I
        title(ax(I,rv), sprintf('%s', drr.labels{rvi}));
        %xlabel(ax(1+jj,rv), sprintf('%s', drr.labels{rvi}));
    end
end

linkaxes(ax);

xlabel(ax(2,1), 'Modulation Frequency (Hz)', 'FontSize', 0.75*fontsize);
ylabel(ax(2,1), 'Frequency (kHz)', 'FontSize', 0.75*fontsize);

% set(ax, 'FontSize', 0.75*fontsize);
% ylabel(ax(2,1), 'Speaker Number', 'FontSize', fontsize)

% Set the colorbar
max_caxis = arrayfun(@(AX) max(caxis(AX)), ax, 'UniformOutput', 1);
max_caxis = max(max_caxis(:));
arrayfun(@(AX) caxis(AX, [0, max_caxis]), ax, 'UniformOutput', 1);

cbar = colorbar(ax(1,end));
cbar.Position = [0.9226 0.6150 0.0146 0.2869];

% ax_ = ax';
aux.abc(ax(:,1), 'location', 'northwestoutside');



%% CF histogram
figh = figure(60+fignum);
clf;

fontsize = 32;

valid_idx = logical( squeeze(prod(~isnan(sum(H(:,1:n_drr,:),1)),2)) );
assert(size(H,3) == nnz(valid_idx), '--> Some of the loaded units are INVALIDE...');


[hist_nodes, hist_edges] = histcounts(1e-3*tbl_valid.CF, 1e-3*f);
hist_edges = 0.5 * (hist_edges(1:end-1) + hist_edges(2:end));
% [] = hist(tbl_valid.CF, 25);
barh = bar(log10(hist_edges), hist_nodes);
if strcmpi('SU' ,finfo.type)
    legend_str = sprintf('%d units', length(tbl_valid.CF));
else
    legend_str = sprintf('%d meas.', length(tbl_valid.CF));
end
legend_h = legend( legend_str, 'Location', 'northwest' );
legend_h.FontSize = 48;
xticks = get(gca, 'XTick');
set(gca, 'XTickLabel', arrayfun(@(X) sprintf('%.1f', X), 10.^xticks, 'UniformOutput', 0));
ylabel('Count');
xlabel('Frequency (kHz)');

set(gca, 'FontSize', fontsize);
title(gca, sprintf('STRF Best Frequencies (%s)', finfo.type));




%% Plot the reconstructed filters
fignum_ = 60+fignum;
figh = figure(fignum_);
clf;

% slc_drr   = drr.dry;
slc_chunk = 1;

obj_loaded = obj_list{1, slc_chunk};
assert(obj_loaded.binwidth == binwidth, '--> You are using different BINWIDTHs!');

if ~isfield(obj_loaded, 'iscausal')
    % backwards cmpatibility 
    obj_loaded.iscausal = true;
end

% Initialize the reconstruction object
obj = reconstruct_c(binwidth,...
    'binwidth', obj_loaded.binwidth,...
    'lags_ms', obj_loaded.lags_ms, ...
    'f', obj_loaded.f, ...
    'lags_ms', obj_loaded.lags_ms,...
    'iscausal', obj_loaded.iscausal, ...
    ...'inv_type', obj_loaded.inv_type, ...
    'remove_bias', obj_loaded.bias_choice ); 
obj.G = obj_list{1, slc_chunk}.G;

% assert(slc_drr == obj_list{slc_drr, slc_chunk}.train_drr);
assert(obj.n_neurons == obj_list{1, slc_chunk}.n_neurons);

    % Plot same units
    %{
    'ToDo: change this!'
    dummy = load('C:\Users\barzelo\Google Drive\codeOnCloud\Reverberation\.data\SU_slc_neurons_from_xls.mat');
    su_slc_neurons_from_xls = dummy.slc_neurons_from_xls;
    dummy = load('C:\Users\barzelo\Google Drive\codeOnCloud\Reverberation\.data\MUA_slc_neurons_from_xls.mat');
    mua_slc_neurons_from_xls = dummy.slc_neurons_from_xls;
    g_filter_list = intersect(mua_slc_neurons_from_xls , su_slc_neurons_from_xls);
    %}
    subxy = [5, 5];
    load('..\.data\MUA_SU_intersect_neurons_from_xls.mat', 'g_filter_list');
    if strcmpi('MUA', finfo.type)
        [];
    elseif strcmpi('SU', finfo.type)
        g_filter_list = 1:length(g_filter_list);  	
    else
        error('-> Unrecognized measurement type (%s)!!',  finfo.type);
    end
    
ax = obj.plot_reconstruction_filters(1e-3*spec_st.f, fignum_, subxy, g_filter_list);
set(ax, 'FontSize', 24);
label_h = aux.abc(ax, 'fontsize', 36, 'label_type', 'numbers');

title(ax(3,1), sprintf('Reconstructed Filters (%s %d units, %s)', finfo.type, obj.n_neurons, drr.labels{slc_drr}));





