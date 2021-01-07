% plot_STRFs_and_BF.m
%
% Plots the histogram of BF-STRFs 

clc;
fignum = 11;

addpath('../');
FigSetup;


%% Load the STRFs data of SUs
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
fn.path = '..\.data\STRFs\';
% fn.name = 'STRF_SU_(25-Aug-2020)_units(103)_bw(5)ms_algo(regression)_fbands(30)_lags(30)ms_cau(1)_trainDRR(3).mat';
fn.name = 'STRF_SU_(27-Aug-2020)_units(103)_bw(1)ms_algo(regression)_fbands(30)_lags(20)ms_cau(1)_trainDRR(3)';
fn.file2load = fullfile(fn.path, fn.name);
load(fn.file2load);


%% Plot best frequency (BF) histogram
figure(fignum);
clf;

fontsize = 34;

bf = 1e-3*tbl_strf.bf;
hist(bf, 20);
% % [N, edges] = histcounts(1e-3*tbl_strf.bf, 20);
% [N, edges] = histcounts(1e-3*tbl_strf.bf, 1e-3*spec_st.f(1:4:end));
% hist(N, edges);
xlabel('Frequency (kHz)');
ylabel('Count');
set(gca, 'FontSize', fontsize);
title(sprintf('STRF Best Frequency of SUs (%d units)', length(bf)));



%% Plot arbitrary STRF
idx_strfs= 101 % options: [21, 101]    % STRF # to plot

binwidth = obj.binwidth; % spec_st.binwidth;    % (ms)
lags_ms  = obj.lags_ms;  %binwidth*n_lags;  	% (ms)
iscausal = 1;
f        = spec_st.f;
jk_flag  = 1;
add_margins= 0;
gausswin_sigma = 0.5;   % smoothing the STRF

fontsize = 34;

% Create a fictitious STRF object (for the plotting)
obj = strf_c(binwidth,...
    'f', f,...
    'iscausal', iscausal, ...
    'lags_ms', lags_ms,...
    'jk_flag', jk_flag); 

obj.strf = tbl_strf.strf{ idx_strfs };

figure(5+fignum);
clf;
obj.plot_strf(add_margins, gausswin_sigma);
colorbar('off');
set(gca, 'FontSize', fontsize); 
colormap jet


% get(gcf, 'Position')
set(gcf, 'Position',[228    64   896   852])



%%
[n_freq, n_lags] = size(tbl_strf.strf{1});
assert(numel(spec_st.f) == n_freq, '--> # of frequency bands are not the same!');

binwidth = spec_st.binwidth;    % (ms)
lags_ms = binwidth*n_lags;  	% (ms)
iscausal = 1;
f        = spec_st.f;
jk_flag  = 1;

% Selected STRFs to plot
slc_strfs = [23, 10, 64, 8];
n_strf = length(slc_strfs);

% Create a fictitious STRF object (for the plotting)
obj = strf_c(binwidth,...
    'f', f,...
    'iscausal', iscausal, ...
    'lags_ms', lags_ms,...
    'jk_flag', jk_flag); 


figure(10+fignum);
clf;

fontsize = 20;
sigma = [];    % 2D Gauss Window
ax = [];
for k = 1:n_strf
    ax(k) = subplot(n_strf,1,n_strf-k+1);
    obj.strf = tbl_strf.strf{ slc_strfs(k) };
    obj.plot_strf([], sigma);
    colorbar('off');
    if k > 1
        set(gca, 'YTickLabel', '');
        set(gca, 'XTickLabel', '');
        xlabel('');
        ylabel('');        
    end
end

set(ax, 'FontSize', fontsize);


%% Load responses
data_type    = 'SU'    % {'SU', MUA'}
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


fprintf('--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);
fprintf('-> data_type: %s\n', data_type);


% Get the valid measurements\columns
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
n_units = nnz(slc.valid_neuron_idx);

% Choosing # of units
H = data.H(:,:,slc.valid_neuron_idx);


%% Show selected responses and their STRF reconstructions
figure(20+fignum);
clf;
fontsize = 20;
addpath(genpath('../fastASD'));

drr     = get_DRR_list_and_indices;
k       = 3;                % select unit # in slc_strfs
unit_idx= slc_strfs(k);  	% STRF's unit # to show
sp      = 6;                % speaker # to plot
fs      = 1/(1e-3*binwidth);


[split_time_idx, n_splits, tbl_metadata] = ... 
    split_spectrogram_into_TIMIT_wav_files(spec_st.binwidth, duration_sec);

t = 1/fs*((split_time_idx(sp,1):split_time_idx(sp,2)-1)' - split_time_idx(sp,1));

obj_array = cell(1, drr.n_drr);
plth      = nan(1, drr.n_drr);
ax        = nan(1, drr.n_drr); 
lgdh      = nan(1, drr.n_drr); 
%scores    = nan(1, drr.n_drr);
for rv = 1:drr.n_drr
    drr_k = drr.ordered(rv);
    X1 = spec_st.Sft{drr_k};  
    y1 = H(:, drr_k, unit_idx);

    [X_train, X_test, y_train, y_test, ~] = train_test_split(X1, y1, ...
        'split_time_idx', split_time_idx,...
        'test_grp', sp ...
    );

    ax(rv) = subplot(drr.n_drr,1,drr.n_drr-rv+1);
    color_k = aux.rpalette(sprintf('new%02d', rv));    
    if strcmpi('SU', data_type)
        bar(t, y_test,...
            'FaceColor', color_k,...
            'EdgeColor', 'none');
        plth(rv) = gca;
    else
        plth(rv) = plot(t, y_test, 'Color', color_k);  % show the response 
    end

    ylabel(ax(rv), drr.labels{drr.ordered(rv)});

    % STRF
    %obj.strf = tbl_strf.strf{ unit_idx };        
    obj_array{rv} = strf_c(binwidth,...
        'f', spec_st.f,...
        'iscausal', 1, ...
        'lags_ms', lags_ms,...
        'algo_type', 'asd',...
        'jk_flag', 0); 
    obj_array{rv}.fit(X_train, y_train);
    y_est = obj_array{rv}.predict(X_test);
    y_est = max(0, y_est);
    gof = goodness(y_test, y_est);
    
    hold on
    plot(t, y_est, 'Color', 'k');
    hold off
    lgdh(rv) = legend(ax(rv), {'PSTH', sprintf('$\\hat{y}_{STRF}$ (CC: %.2f)', gof.CC)});
end
set(ax, 'YTickLabel', '');
set(ax(2:end), 'XTickLabel', '');
xlabel(ax(1), 'Time (sec)');
set(ax, 'FontSize', fontsize);
title(ax(end), sprintf('$BF_{strf}$: %.2f kHz', 1e-3*obj_array{1}.bf));
set(lgdh, 'FontSize', 16);
aux.abc(ax, 'FontSize', 32);




%%
figure(25+fignum);
clf;
fontsize = 20;
ax2 = nan(1, drr.n_drr); 

for k = 1:drr.n_drr
    ax2(k) = subplot(drr.n_drr,1,drr.n_drr-k+1);
    obj_array{k}.plot_strf;
    
end


