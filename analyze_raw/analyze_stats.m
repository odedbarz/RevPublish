%
% analyze_stats.m
%
% Description:
% Calculate various statistics and plot them.
%
% Reference:
%     RMD: response modulation depth
%     MG : modulation gain
%     CCr: correlation coefficient between DRY and other DRR conditions


clc
fignum = 10;
verbose = 1;
setup_environment('../');




%% Plot properties
fontsize = 32;
fontsize_big = 42;
fontsize_bigger = 64;
markersize = 24;




%% Load data
% 
%   Name                Size                 Bytes  Class     Attributes
% 
%   H_units          7200x5x10             2880000  double              
%   fn                  1x1                   2206  struct              
%   obj_list            5x12              13504080  cell                
%   sorted_list       241x1                   1928  double              
%   spec_st             1x1               15039029  struct              
%   splits              1x1                  58552  struct              
%   tbl_data          241x20                339094  table               
%
% Notes:
% >> reconstruct_XXX_(14-Jan-2021)_...: in which XXX is SU or MUA; SUs & MUAs are 
%                                      taken from different recording sites.
% >> reconstruct_XXX_(02-Jun-2021)_...: SUs & MUAs are taken from THE SAME recording sites.
 %
data_type   = 'SU';       % {'SU', MUA'}
fn_path= '../_data/Reconstruct/';
data_type   = upper(data_type);
switch data_type
    case 'SU'
        fn_template = 'reconstruct_SU_(07-Jun-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';       
        unit_list = 100;    % load max number of available units
        
    case 'MUA'
        fn_template = 'reconstruct_MUA_(07-Jun-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';        
        unit_list = 100;    % load max number of available units

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end

% Loading data & metadata
fn_n = sprintf(fn_template, unit_list(1));
aux.cprintf('String', '\n-> Loading FIRST file to get preliminary data <%s>...\n', fn_n);
warning off
data = load( fullfile(fn_path, fn_n) );
warning on



%% Initialization
drr = get_DRR_list_and_indices; 
n_drr = drr.n_drr;
drr_idx = drr.ordered;
% len_unit_list = length(unit_list);

splits   = data.splits;
spec_st  = data.spec_st;
stim_st  = data.stim_st;
H_sorted = data.H_sorted;           % size(H_sorted) = [time\samples x DRR x units]
sorted_list = data.sorted_list;     % the order in which the neurons are sorted in H_sorted
tbl_data = data.tbl_data;
n_splits = splits.n_splits;
n_bands  = spec_st.n_bands;
n_units  = size(H_sorted,3);
assert( unit_list == n_units, 'ERROR: there is a missfit between DESIRED and LOADED number of units');
 

% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)



%% UNSORT order in H_sorted
% * NOTE: H_sorted is *NOT* arranged by drr.ordered (dim 2) 
[~, idx_reset_sort] = sort( sorted_list );
H = H_sorted(:, drr_idx, idx_reset_sort);




%% BROAD-BAND CC_broadband_envelope
% Calculate the broadband CCs
% ! unused
%{
win_env_lpf_ms = 60;
fs_dwnsmp      = fs_timeaxis;
[R, args]      = broadband_envelopes_CCs(stim_st.Y, stim_st.fs, win_env_lpf_ms, fs_dwnsmp);
CC_broadband_envelope = R(drr.dry, drr.sortby(1:n_drr));
%}




%% Get Best-Frequencies (BF)
% Loads the BF table:
%
% tbl_BF:
% 
%     neuron      BF        BF_cc       BF_pv  
%     ______    ______    ________    ________
%        1        5898     0.21813     0.21813
%        2      5325.8     0.64889     0.64889
%       ...       ...        ...         ...
%    
fn_path_data = '../_data';

switch data_type
    case 'SU'
        fn_data = 'data_SU_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';  
        fn_strf = 'STRF_SU_(22-Apr-2021)_units(103)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(1)_trainDRR(3).mat';

    case 'MUA'
        fn_data = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';  
        fn_strf = 'STRF_MUA_(22-Apr-2021)_units(241)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(1)_trainDRR(3).mat';

    otherwise
        error('--> Unrecognized DATA_TYPE!');
end

% [~, tbl_BF] = find_best_unit_set('CC', 'fn', fullfile(fn_path_data, fn_data) );
% BF = tbl_BF.BF;
    
% Get spike-units for both SU & MUA
file_path = load.path_to_data('_data');
file_template = 'data_%s_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
sorted_list = find_best_unit_set('SPK', 'fn_template', ...
    {file_path, file_template, data_type});    % {'SPK', 'NOSPK'}

% Load the STRF table with the BFs
file_path = load.path_to_data('Analysis');
data = load( fullfile(file_path, fn_strf) );

% Get the spike-units's BFs
BF = data.tbl_strf.bf(sorted_list);




%% Correlation Coefficients between STIMULI
Sdry = spec_st.Sft{drr.dry};
CCs = nan(1, n_drr);    % CCs of responses
for k = 1:n_drr
    rv = drr.ordered(k);
    Sk = spec_st.Sft{rv};
    
    % CCs: correlation between SRY & DRR stimuli    
    gof = goodness(Sdry, Sk);
    CCs(k) = gof.CC;
end




%% Plot: RMD, MG, Kurtosis, and CCr
% RMD: response modulation depth
% MG : modulation gain
% CCr: correlation coefficient between DRY and other DRR conditions

% Select data for analysis
scores.CCr  = nan(n_units, n_drr);       
scores.MDr  = nan(n_units, n_drr);
scores.MG   = nan(n_units, n_drr);
scores.skew = nan(n_units, n_drr);
scores.ku   = nan(n_units, n_drr);

% The DRY response 
Hdry = squeeze( H(:, 1, :) );        % size(H_sorted) = [time\samples x DRR x units]

for n = 1:n_units
    %In = sorted_list(n);
    In = n;
    
    % Stimulus envelope; find the closest frequency band to the neuron's CF 
    [~, idx_cf] = min(abs(f - BF(In)));

    for k = 1:n_drr                
         % CC
        CC = corrcoef(Hdry(:,n), H(:, k, n));
        scores.CCr(n,k) = CC(1,2);

        % Response MD (modulation depth)
        scores.RMD(n,k) = MD( H(:, k, n), 1 );

        % Envelope's MD 
        rv = drr.ordered(k);
        env_db = spec_st.Sft{rv}(idx_cf,:);        
        %cfloor = spec_st.db_floor;              % The spectrogram (Sft) amplitudes 
        %cmax = spec_st.max_Sft;                 % are of log scale, so we need to 
        %env_ = cmax*10.^((env_db + cfloor)/20); % revert it
        %MD_env = MD( env_ );
        MD_env = MD( env_db );

        % MG (modulation gain) = RMD/MDenv
        scores.MG(n,k) = db( scores.RMD(n,k)./MD_env );
        scores.skew(n,k) = skewness( H(:, k, n) );
        scores.kr(n,k) = kurtosis( H(:, k, n) );
    end
     
end




%% RMD
figure(100+fignum);
clf;
ax = gca;
plot_dotbox(scores.RMD, 'labels', drr.labels(drr_idx));
set(ax, 'FontSize', fontsize);
ylabel(aux.ctitle('RMD', '$(\sqrt{2}\sigma/\mu)$'));
xlabel('Direct to Reverberation Ratio (dB)');
title(sprintf('%d %s', n_units, data_type));

% Wilcoxon signed rank test between SU & MUA 
tbl_pv = array2table(nan(2, n_drr-1), 'VariableNames',...
    {'dryTOdB9_4', 'dB9_4TOdB4_8', 'dB4_8TOdBm2_5', 'dBm2_5TOdBm8_2'} );
for k = 1:n_drr-1
    [tbl_pv{1,k}, tbl_pv{2,k}] = signrank(scores.RMD(:,k), scores.RMD(:,k+1));
end
fprintf('\n RMD: Wilcoxon signed rank test between SU (reverberant vs. dry speech)\n')
tbl_pv




%% MG
figure(105+fignum);
clf;
ax = gca;
plot_dotbox(scores.MG, 'labels', drr.labels(drr_idx));
set(ax, 'FontSize', fontsize);
ylabel('MG (dB)');
xlabel('Direct to Reverberation Ratio (dB)');
title(sprintf('%d %s', n_units, data_type));

% Wilcoxon signed rank test between SU & MUA 
tbl_pv = array2table(nan(2, n_drr-1), 'VariableNames',...
    {'dryTOdB9_4', 'dB9_4TOdB4_8', 'dB4_8TOdBm2_5', 'dBm2_5TOdBm8_2'} );
for k = 1:n_drr-1
    [tbl_pv{1,k}, tbl_pv{2,k}] = signrank(scores.MG(:,k), scores.MG(:,k+1));
end
fprintf('\n MG: Wilcoxon signed rank test between SU (reverberant vs. dry speech)\n')
tbl_pv




%% CC
figure(110+fignum);
clf;
ax = gca;
plot_dotbox(scores.CCr, 'labels', drr.labels(drr_idx));
set(ax, 'FontSize', fontsize);
ylabel('CC');
hold on
plth = plot(CCs, 'sk:', 'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k');
hold off
legend(plth, 'CCs stimulus');
xlabel('Direct to Reverberation Ratio (dB)');


% Wilcoxon signed rank test between SU & MUA 
tbl_pv = array2table(nan(2, n_drr-1), 'VariableNames',...
    {'dryTOdB9_4', 'dB9_4TOdB4_8', 'dB4_8TOdBm2_5', 'dBm2_5TOdBm8_2'} );
for k = 1:n_drr-1
    [tbl_pv{1,k}, tbl_pv{2,k}] = signrank(scores.CCr(:,k), scores.CCr(:,k+1));
end
fprintf('\n CC: Wilcoxon signed rank test between SU (reverberant vs. dry speech)\n')
tbl_pv






%%
figure(120+fignum);
clf;
ax = gca;
plot_dotbox(scores.kr, 'labels', drr.labels(drr_idx));
set(ax, 'FontSize', fontsize);
ylabel('Kurtosis');
ylim([min(ylim), 15]);
xlabel('Direct to Reverberation Ratio (dB)');
ylim([1.5, 10]);
title(sprintf('%d %s', n_units, data_type));



%%
figure(125+fignum);
clf;
ax = gca;
plot_dotbox(scores.skew, 'labels', drr.labels(drr_idx));
set(ax, 'FontSize', fontsize);
ylabel('Kurtosis');
ylim([min(ylim), 15]);
xlabel('Direct to Reverberation Ratio (dB)');
% ylim([1.5, 10]);
title(sprintf('%d %s', n_units, data_type));




























