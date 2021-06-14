% analyze_CC_RMD_MG.m


clc
fignum = 11;
verbose = 1;

setup_environment('../../');
  
drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);




%% Load a list of sorted units 
% (sorted by MUA repitability between raw measurements)
fn_unit_list = fullfile( load.path_to_data('Analysis'), 'idx_MUA_good_sorted_unit_thr(0-7).mat');
dummy = load(fn_unit_list);
sorted_unit_list = dummy.unit_list;



%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type = 'MUA';       % {'SU', MUA'}
data_type = upper(data_type);

switch data_type
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        fn.load.file = 'data_SU_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        unit_list = 103; %[10, 25, 50, 103];
        
        %fn.strf.path = load.path_to_data('Analysis');
        %fn.strf.file = 'STRF_SU_(22-Apr-2021)_units(103)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(1)_trainDRR(3).mat';
        

    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        fn.load.file = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';  
        unit_list = 241; %[10, 25, 50, 103, 241];
        
        %fn.strf.path = load.path_to_data('Analysis');
        %fn.strf.file = 'STRF_MUA_(22-Apr-2021)_units(241)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(1)_trainDRR(3).mat';

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end

% Load the data
fn.load.path = load.path_to_data('data');
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data = load(fn.load.fullfile);

% Load tbl_strf
% dummy = load(fullfile(fn.strf.path, fn.strf.file));
% tbl_strf = dummy.tbl_strf;
[~, tbl_BF] = find_best_unit_set('CC', 'fn', fn.load.fullfile);  % {'CC'}


spec_st   = data.spec_st;
stim_st   = data.stim_st;
tbl_data  = data.(sprintf('tbl_%s', data_type));

n_units   = height(tbl_data);     % total available units
len_unit_list = length(unit_list);
duration_sec = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);




%% Reconstruction parameters
iscausal       = 0;                    % use reconstructed causal filters?
lags_ms        = 30;                   % (ms) maximum lags
binwidth       = spec_st.binwidth;     % (ms)
n_bands        = spec_st.n_bands;
win_size_ms    = spec_st.win_size_ms; 
jackknife_flag = 1;     
algo_type      = 'regression';    % {'regression', 'asd', 'svd'}
duration_sec   = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

if verbose
    aux.cprintf('UnterminatedStrings', '\n    Data:\n');
    aux.cprintf('UnterminatedStrings', '--> data_type   : %s\n', data_type);
    aux.cprintf('UnterminatedStrings', '--> unit_list   : [%s]\n', num2str(unit_list));
    aux.cprintf('UnterminatedStrings', '--> n_units     : %g (all available units)\n', n_units);
    aux.cprintf('UnterminatedStrings', '--> duration_sec: %g ms\n', duration_sec);
    aux.cprintf('UnterminatedStrings', '    Reconstruction:\n');
    aux.cprintf('UnterminatedStrings', '--> causality   : %d\n', iscausal);
    aux.cprintf('UnterminatedStrings', '--> lags_ms     : %g ms\n', lags_ms);
    aux.cprintf('UnterminatedStrings', '--> binwidth    : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands     : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms : %g ms\n', win_size_ms);
    aux.cprintf('UnterminatedStrings', '--> is_jackknife: %d\n', jackknife_flag);
end




%% Correlation Coefficients between STIMULI
CCs = CC_stimuli( spec_st );



%% Analyze for RMD, MG, CCer, and KURTOSIS
% RMD: response modulation depth
% MG : modulation gain
% CCer: correlation coefficient between DRY envelope (on BF) and resposnes
n_units = size(data.H,3);

% Select data for analysis
H = squeeze( data.H(:,:, 1:n_units) );

clear scores
scores.CCer = nan(n_units, n_drr);       % % CC(response dry & response)
scores.MDr  = nan(n_units, n_drr);
scores.MG   = nan(n_units, n_drr);
scores.RMD  = nan(n_units, n_drr);
scores.skew = nan(n_units, n_drr);
scores.ku   = nan(n_units, n_drr);


for n = 1:n_units
    %un = sorted_unit_list(n);
    
    % Option #1:
    % Stimulus envelope; find the closest frequency band to the neuron's CF 
    [~, idx_bf] = min(abs(spec_st.f - tbl_BF.BF(n)));
    
    % Option #2
    % Stimulus envelope: best correlation between envelope and response
    %Sft_ = zca(spec_st.Sft{drr.dry}, 2);
    %[~, idx_bf] = max( Sft_ * H(:,drr.dry,n) );
    %idx_bf = 30-idx_bf+1;
    
    % Option #3:
    %[~, idx_bf] = min(abs(spec_st.f - tbl_data.CF(n)));        
    
    % Option #3:    % '####### DEBUG ######'
    % Random envelope
    %idx_bf = randi(length(f));
    

    % Signal Envelope
    yenv = spec_st.Sft{ drr.dry }(idx_bf,:); 
    
    for k = 1:n_drr
        rv = drr.ordered(k);
        
        % Response MD (modulation depth)
        scores.RMD(n,k) = MD(H(:,rv,n), 1);
        
        % Envelope's MD 
        % Option #1: using the spectrogram
        %yenv = spec_st.Sft{rv}(idx_bf,:);         
        %{
        cfloor = spec_st.db_floor;              % The spectrogram (Sft) amplitudes 
        cmax = spec_st.max_Sft;                 % are of log scale, so we need to 
        yenv_10 = cmax*10.^((yenv + cfloor)/20); % revert it
        MD_env = MD( yenv_10 );
        %}
        %
        % Option #2: using the stimulus (more computation time)
        % %{
        %stim_env_rv = calc_stimulus_envelope(stim_st.Y(:,rv), tbl_valid.CF(n), fs, stim_st.fs);
        stim_env_rv = calc_stimulus_envelope(stim_st.Y(:,rv), spec_st.f(idx_bf), ...
            1/(1e-3*binwidth), stim_st.fs);
        MD_env = MD( stim_env_rv, 1 );
        %}
        
        % MG (modulation gain) = RMD/MDenv
        scores.MG(n,k) = db(scores.RMD(n,k)./MD_env);
        
        % Kurtosis & Skewness
        scores.skew(n,k) = skewness( H(:,rv,n) );
        scores.kr(n,k) = kurtosis( H(:,rv,n) );

        % CCer: DRY(response)-to-DRR(response)
        %dummy = corrcoef(H(:,drr.dry,n), H(:,rv,n));
        if strcmpi(data_type, 'SU')
            su_win_size = 1;   % su_win_size * binwidth; su_win_size == 1 means NO window
            su_win = hann(su_win_size)/norm(hann(su_win_size));
            su_win = su_win(ceil(end/2):end);
            hn = filter(su_win, 1, H(:,rv,n));
            %hn = H(:,rv,n);
        else
            hn = H(:,rv,n);
        end
        dummy = corrcoef(yenv, hn);
        scores.CCer(n,k) = dummy(1,2);
        
     end
end




%%
% [~, cc_sorted_idx] = sort(scores.CCer(:,1), 'descend');
switch upper(data_type)
    case 'SU'
        cc_sorted_idx = 1:n_units;      
        
    case 'MUA'
        cc_sorted_idx = sorted_unit_list;
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end

scores.CCer_sorted = scores.CCer(cc_sorted_idx, :);
scores.RMD_sorted  = scores.RMD(cc_sorted_idx, :);
scores.MG_sorted   = scores.MG(cc_sorted_idx, :);



%% Load shared measurements (between SU & MUA)
% file_template = 'data_%s_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
% 
% sorted_list = find_best_unit_set('SPK', 'fn_template', ...
%     {load.path_to_data('_data'), file_template, data_type});    % {'SPK', 'NOSPK'}
% 
% 
% CCer_sorted = scores.CCer_sorted(sorted_list,:);
% n_sorted_units_to_plot = length(sorted_list);

% *** OR ***
n_sorted_units_to_plot = 100;
CCer_sorted = scores.CCer_sorted(1:n_sorted_units_to_plot,:);


%% === SORTed ===
% CCer_SORTed: correlation coefficients between bandpass-envelope and response
figure(0+fignum);
% clf;

markersize = 32;
fontsize = 28;
fontsize_big = 42;
fontsize_bigger = 48;

switch upper(data_type)
    case 'SU'
        ax = subplot(1,2,1);        
        
    case 'MUA'
        ax = subplot(1,2,2);
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end

% Add a violine plot
% h = aux.violinplot(CCer_sorted, {}, 'ShowMean', true, 'ShowData', false);
h = aux.violinplot(CCer_sorted, drr.labels(drr.ordered), 'ShowMean', true, 'ViolinAlpha', 0.5); %, 'ShowData', false);
for k = 1:n_drr
    h(k).ViolinColor = aux.rpalette(k);
end

plth = [];
% pars = plot_dotbox(CCer_sorted, 'labels', drr.labels(drr.ordered));
% plth = pars.points_h;
set(ax, 'FontSize', fontsize);
hold on
plth(end+1) = plot(CCs, 'sk:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k');
hold off
% legend(plth, {'CC (Broadband Envelope)', 'CC (Stimulus DRY-to-DRR)'});
% legend(plth(1:end-1), {'$S_{dry}$ vs. $\hat{S}_{drr}$', '$S_{drr}$ vs. $\hat{S}_{drr}$'}, 'FontSize', fontsize_big);
% aux.abc(ax, 'fontsize', 2*fontsize, 'location', 'northwestoutside');

switch upper(data_type)
    case 'SU'
        ylabel('CC (Envelope-to-Response)', 'FontSize', fontsize_big);
        set(gca, 'Position', [0.0970, 0.1100, 0.3928, 0.8038]);

    case 'MUA'
        set(ax, 'YTickLabel', '');
        set(gca, 'Position', [0.5279, 0.1100, 0.3928, 0.8038]);

end
set(gca, 'Box', 'On');
xlabel('DRR', 'FontSize', fontsize_big);
title(sprintf('%d %ss', n_sorted_units_to_plot, data_type), ...
    'FontSize', fontsize_bigger);
ylim([-0.15, 1.1]);


%% Plot RMD\MG\CC\Kurtosis
% %{
figure(10+fignum);
clf;

fontsize = 38;
markersize = 32;

ax = subplot(2,1,1);
% plot_dotbox(scores.RMD, 'labels', drr.labels(drr.ordered));
plot_dotbox(scores.RMD_sorted(1:n_sorted_units_to_plot,:), 'labels', drr.labels(drr.ordered));
ylabel(aux.ctitle('RMD', '$(\sqrt{2}\sigma/\mu)$'));
set(gca, 'XTickLabel', '');
xlabel('');
title(ax(1), sprintf('%s', data_type));

ax(2) = subplot(2,1,2);
% plot_dotbox(scores.MG, 'labels', drr.labels(drr.ordered));
plot_dotbox(scores.MG_sorted(1:n_sorted_units_to_plot,:), 'labels', drr.labels(drr.ordered));
ylabel('MG (dB)');
set(gca, 'XTickLabel', '');
xlabel('');

% ax(3) = subplot(3,1,3);
% plot_dotbox(scores.CCer_sorted(1:n_sorted_units_to_plot,:), 'labels', drr.labels(drr.ordered));
% ylabel('$CC$');

% hold on
% plth = plot(CCs, 'sk:', 'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k');
% hold off
% legend(plth, 'CCs stimulus');
set(ax, 'FontSize', fontsize);
xlabel('Direct to Reverberation Ratio (dB)', 'FontSize', fontsize);
% aux.abc(ax, 'fontsize', 2*fontsize, 'location', 'northwestoutside');
%}




%% errorbar RMD\MG\CC\Kurtosis
% %{
figure(20+fignum);
clf;

fontsize = 38;
markersize = 32;

ax = subplot(2,1,1);
boxplot(scores.RMD_sorted(1:n_sorted_units_to_plot,:), 'labels', drr.labels(drr.ordered))
ylabel(aux.ctitle('RMD', '$(\sqrt{2}\sigma/\mu)$'));
set(gca, 'XTickLabel', '');
xlabel('');
title(ax(1), sprintf('%d %s', n_sorted_units_to_plot, data_type));

ax(2) = subplot(2,1,2);
boxplot(scores.MG_sorted(1:n_sorted_units_to_plot,:), 'labels', drr.labels(drr.ordered))
ylabel('MG (dB)');
set(gca, 'XTickLabel', '');
xlabel('');

% ax(3) = subplot(3,1,3);
% boxplot(scores.CCer_sorted(1:n_sorted_units_to_plot,:), 'labels', drr.labels(drr.ordered))
% ylabel('$CC$');

% hold on
% plth = plot(CCs, 'sk:', 'MarkerSize', 0.4*markersize, 'MarkerFaceColor', 'k');
% hold off
legend(plth, 'CCs stimulus');
xlabel('Direct to Reverberation Ratio (dB)');
set(ax, 'FontSize', fontsize);
aux.abc(ax, 'fontsize', 2*fontsize, 'location', 'northwestoutside');
% ylim([min(ylim), 1]);
%}







