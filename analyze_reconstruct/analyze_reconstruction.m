%
% analyze_reconstruction.m
%
% Description:
% Compares spectrogram reconstructions. The comparison is between DRY and other
% DRR conditions.
%
%

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
data_type   = 'MUA';       % {'SU', MUA'}
fn_path = load.path_to_data('Reconstruct');

data_type   = upper(data_type);
switch data_type
    case 'SU'
        fn_template = 'reconstruct_SU_(14-Jan-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';       
        unit_list = [10, 25, 50, 103];
        
    case 'MUA'
        fn_template = 'reconstruct_MUA_(14-Jan-2021)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3).mat';        
        unit_list = [10, 25, 50, 103, 150, 241];

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end




%% Initialization
drr = get_DRR_list_and_indices; 
n_drr = drr.n_drr;
drr_idx = drr.sortby(1:n_drr);
len_unit_list = length(unit_list);


% Loading first file to get preliminary metadata
%        H_sorted: [7200×5×10 double]
%              fn: [1×1 struct]
%        obj_list: {5×12 cell}
%     sorted_list: [1×241 double]
%         spec_st: [1×1 struct]
%          splits: [1×1 struct]
%         stim_st: [1×1 struct]
%        tbl_data: [241×20 table]
fn_n = sprintf(fn_template, unit_list(1));
aux.cprintf('String', '\n-> Loading FIRST file to get preliminary data <%s>...\n', fn_n);
warning off
dummy = load( fullfile(fn_path, fn_n), 'splits', 'spec_st', 'stim_st', 'tbl_data' );
warning on

splits   = dummy.splits;
spec_st  = dummy.spec_st;
stim_st  = dummy.stim_st;
tbl_data = dummy.tbl_data;
n_splits = splits.n_splits;
n_bands  = spec_st.n_bands;
n_units  = height(tbl_data);

% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs_timeaxis = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)




clear scores  %sti sti2
scores.n_drr        = n_drr;
scores.info         = 'X_dry vs. X_est';
scores.n_units      = nan(1, len_unit_list);
scores.CC           = nan(n_drr, n_splits, len_unit_list);
scores.mse          = nan(n_drr, n_splits, len_unit_list);
scores.nmse         = nan(n_drr, n_splits, len_unit_list);
scores.CCf          = nan(n_drr, n_splits, len_unit_list, n_bands);

% For comparing between stimuli of various DRR conditions and the
% reconstructed spectrograms
scores2           = scores;
scores2.info      = 'X_drr vs. X_est';




%%
for n = 1:len_unit_list
    % Load the file
    fn_n = sprintf(fn_template, unit_list(n));
    aux.cprintf('Comments', '\n-> Loading <%s>...\n', fn_n);
    warning off
    data = load(fullfile(fn_path, fn_n));
    warning on
    
    obj_list= data.obj_list;
    spec_st = data.spec_st;
    splits = data.splits;
    
    % # of neurons
    n_neurons = obj_list{1}.n_neurons;
    assert( unit_list(n) == n_neurons, 'ERROR: something is wrong!' );

    if verbose
        fprintf('--> n_neurons  : %d\n', n_neurons);
        fprintf('--> binwidth   : %d\n', binwidth);
        fprintf('--> n_bands    : %d\n', n_bands);
        fprintf('--> win_size_ms: %d\n', win_size_ms);
        fprintf('--> n_drr      : %d\n', n_drr);
        fprintf('--> n_splits   : %d\n', n_splits);
    end
    

    % Loop over the splits (exclusive intervals)
    scores.n_units(n) = n_neurons;
    for sp = 1:n_splits 
        % Cut out the testing chunk of the spectrogram
        X_dry = spec_st.Sft{drr.dry}(:, sp == splits.idx);    
        
        for k = 1:n_drr    
            rv = drr.ordered(k);
            
            % Cut out the testing chunk of the spectrogram
            X_est                = obj_list{rv,sp}.X_est;
            gof                  = goodness(X_dry, X_est);    
            scores.CC(k,sp,n)    = gof.CC;
            scores.mse(k,sp,n)   = gof.mse;
            scores.nmse(k,sp,n)  = gof.nmse;
            scores.CCf(k,sp,n,:) = gof.CCf;            
            
            % For comparing between stimuli of various DRR conditions and the
            % reconstructed spectrograms
            X_drr                 = spec_st.Sft{rv}(:, sp == splits.idx);            
            gof2                  = goodness(X_drr, X_est);  
            scores2.CC(k,sp,n)   = gof2.CC;
            scores2.mse(k,sp,n)  = gof2.mse;
            scores2.nmse(k,sp,n) = gof2.nmse;
            scores2.CCf(k,sp,n,:)= gof2.CCf;
        end
    end
    

end
 

% Statistics
%
% Dims: [drr x # splits x # units]
%
% DRY vs RECONSTRUCTION
% idx_dry_ordered = 1;    % all measurements are already ORDERED by DRR.ORDERED
scores.mu.CC    = squeeze( median(scores.CC, 2) );
scores.SE.CC    = squeeze( mad(scores.CC,[],2)/sqrt(n_splits) );
scores.mu.mse   = squeeze( median(scores.mse, 2) );
scores.SE.mse   = squeeze( mad(scores.mse,[],2)/sqrt(n_splits) );
scores.mu.nmse  = squeeze( median(scores.nmse, 2) );
scores.SE.nmse  = squeeze( mad(scores.nmse,[],2)/sqrt(n_splits) );
scores.mu.CCf   = squeeze( median(scores.CCf(:,:,len_unit_list,:), 2) );
scores.SE.CCf   = squeeze( mad(scores.CCf(:,:,len_unit_list,:), [], 2 ) );

% DRR vs RECONSTRUCTION
scores2.mu.CC   = squeeze( median(scores2.CC, 2) );
scores2.SE.CC   = squeeze( mad(scores2.CC,[],2)/sqrt(n_splits) );
scores2.mu.mse  = squeeze( median(scores2.mse, 2) );
scores2.SE.mse  = squeeze( mad(scores2.mse,[],2)/sqrt(n_splits) );
scores2.mu.nmse = squeeze( median(scores2.nmse, 2) );
scores2.SE.nmse = squeeze( mad(scores2.nmse,[],2)/sqrt(n_splits) );
scores2.mu.CCf   = squeeze( median(scores2.CCf(:,:,len_unit_list,:), 2) );
scores2.SE.CCf   = squeeze( mad(scores2.CCf(:,:,len_unit_list,:), [], 2 ) );



%% BROAD-BAND CC_broadband_envelope
% Calculate the broadband CCs
% ! unused
%{
win_env_lpf_ms = 60;
fs_dwnsmp      = fs_timeaxis;
[R, args]      = broadband_envelopes_CCs(stim_st.Y, stim_st.fs, win_env_lpf_ms, fs_dwnsmp);
CC_broadband_envelope = R(drr.dry, drr.sortby(1:n_drr));
%}


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




%% CC vs. various DRRs (bar plot)
figure(0+fignum);
clf;
    

x = 1:n_drr;
M = [scores.mu.CC(:,end), scores2.mu.CC(:,len_unit_list)];
errbar = [scores.SE.CC(:,len_unit_list), scores2.SE.CC(:,len_unit_list)];
h = bar(M); 
set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr_idx));
h(1).FaceColor = aux.rpalette(sprintf('new%02d',1));
h(2).FaceColor = aux.rpalette(sprintf('new%02d',2));

dx = 0.2*h(1).BarWidth;
hold on
% ERROR BARs
plth = errorbar([(1:5)'-dx, (1:5)'+dx], M, errbar);

% CCs
plot(1:n_drr, CCs, 'ks:',...
    'MarkerSize', markersize, 'MarkerFaceColor', 'k' );

set(gca, 'FontSize', fontsize);

hold off
for k = 1:length(plth)
    plth(k).Color = [0 0 0];                            
    plth(k).LineStyle = 'none'; 
    plth(k).LineWidth = 2;
end

ylim([0.0, 1.4]);
legend(h, {'Compare with Dry Condition', 'Compared with Same DRR condition'}, 'FontSize', fontsize_big);
ylabel('CC', 'FontSize', fontsize_big);
xlabel('DRR', 'FontSize', fontsize_big);
title( sprintf('%d %ss', n_units, data_type) );
aux.ctitle('Dry Stimulus vs. Reconstructions', sprintf('(%d %ss)', n_units, data_type));





%% CC vs. Speaker sex 
figure(5+fignum);
clf;
ax = gca;


D = scores.CC(:,:,len_unit_list)';
plth = plot(D, 's:', 'MarkerSize', markersize);
if strcmpi('MUA', data_type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end
set(gca, 'FontSize', fontsize);

xlabel('Speaker Sex', 'FontSize', fontsize_big);
ylabel('CC', 'FontSize', fontsize_big);
aux.ctitle('Dry Stimulus vs. Reconstructions', sprintf('(%d %ss)', n_units, data_type));

legend( drr.labels{drr_idx}, 'Location', 'southeast', 'FontSize', fontsize_big);
axis tight 
ylim([0, 1]);

speaker_sex_list = {'M$_1$','M$_2$','M$_3$','F$_1$','M$_4$','F$_2$','F$_3$',...
    'F$_4$','M$_5$','F$_5$','F$_6$','M$_6$'};
xticks = get(ax(1), 'XTick');
set(ax(1), 'XTickLabel', speaker_sex_list(xticks));




%% Correlations between frequencies of selected intervals (jackknife)
figure(10+fignum);
clf;

ax = gca;

plth = errorbar(1e-3*f, scores.mu.CCf, scores.SE.CCf, 's-', 'MarkerSize', markersize);
if strcmpi('MUA', data_type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end
ylim([0.75, 1.0]);

set(gca, 'FontSize', fontsize);
xlabel('Frequency (kHz)', 'FontSize', fontsize_big);
ylabel('CC', 'FontSize', fontsize_big);
aux.ctitle('Dry Stimulus vs. Reconstructions', sprintf('(%d %ss)', n_units, data_type));





%% Spectrograms: Stimuli vs. Reconstructions
figh = figure(15+fignum);
clf;

spec_thr    = 20;       % (dB) threshold the spectrograms from below for clarity 
color_axis  = [40, 90];  % (dB) set same range of colors to all spectrograms 

speaker = 6;      % split number \ jackknife number to use
nt = nnz(speaker == splits.idx);
t = spec_st.t(1:nt);       % (sec)
ax = zeros(n_drr,2);

for k = 1:n_drr   
    rv = drr.ordered(k);
    
    % Stimulus
    Sft = spec_st.Sft{rv}(:, speaker == splits.idx);   
    Sft_thr = max(spec_thr, Sft);    
    
    % Reconstruction 
    Sft_est  = obj_list{rv,speaker}.X_est;
    Sft_est_thr = max(spec_thr, Sft_est);
    
    %ax(k,1) = subplot(5,2,1+2*(k-1));            
    ax(k,1) = subplot(5,2, 5*2 - (1+2*(k-1)));    
    spec.plot_spectrogram(t, 1e-3*f, Sft_thr, figh, 1);
    
    %ax(k,2) = subplot(5,2,2*k);
    ax(k,2) = subplot(5,2, 5*2-2*(k-1));
    spec.plot_spectrogram(t, 1e-3*f, Sft_est_thr, figh, 1);
    
    % Add the DRRs to the ylabels
    if 3 ~= k
        ylabel(ax(k,1), drr.labels{rv});
    else
        ystr = aux.ctitle('Frequency (kHz)\\', drr.labels{rv});
        ylabel(ax(k,1), ystr);        
    end

end

ax = flipud(ax);
arrayfun(@(X) caxis(X, color_axis), ax);   % (dB) set the range of colors

set(ax, 'FontSize', fontsize);
set(ax(1:4,:), 'XTickLabel', '');
set(ax(:,2), 'YTickLabel', '');
xlabel(ax(end,1), 'Time (sec)', 'FontSize', fontsize_big);
xlabel(ax(end,2), 'Time (sec)', 'FontSize', fontsize_big);
title(ax(1,1), 'Stimulus');
title(ax(1,2), sprintf('Reconstruction (%d %ss)', n_units, data_type));
aux.abc(ax(1,:), 'fontsize', fontsize_bigger, 'location', 'northwestoutside');
linkaxes(ax);       
            
% Add a colorbar and make sure that the axes doesn't change in size (push 
% the colorbar "outside").
if 1
    ax_j = ax(end,end);
    pos = get(ax_j, 'Position');
    
    colorbar(ax_j);
    set(ax_j, 'Position', pos);
    
end































