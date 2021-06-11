%
% main_BFs.m
%


clc
fignum = 11;
verbose = 1;

setup_environment('../../');
  
drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);




%% Load data
%   Run [main_STRF_loopover.m] again to update this file if needed
% 
%   Name            Size               Bytes  Class     Attributes
% 
%   fn              1x1                 2170  struct              
%   spec_st         1x1             15039029  struct              
%   splits          1x1                58376  struct              
%   strf_st         1x1                 7365  struct              
%   tbl_strf      241x9             75466618  table               
data_type   = 'SU';       % {'SU', MUA'}
data_type   = upper(data_type);

switch data_type
    case 'SU'
        fn.load.file = 'STRF_SU_(22-Apr-2021)_units(103)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(1)_trainDRR(3).mat';       

    case 'MUA'
        fn.load.file = 'STRF_MUA_(22-Apr-2021)_units(241)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(1)_trainDRR(3).mat';

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end

fn.load.path = load.path_to_data('Analysis');
fn.load.fullfile = fullfile( load.path_to_data('Analysis'), fn.load.file );
load(fn.load.fullfile);


n_units   = height(tbl_strf);     % total available units
duration_sec = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);




%% Reconstruction parameters
iscausal       = strf_st.iscausal;      % use reconstructed causal filters?
lags_ms        = strf_st.lags_ms;       % (ms) maximum lags
binwidth       = spec_st.binwidth;      % (ms)
n_bands        = strf_st.n_bands;       % (ms)

if verbose
    aux.cprintf('UnterminatedStrings', '\n    Data:\n');
    aux.cprintf('UnterminatedStrings', '--> data_type   : %s\n', data_type);
    aux.cprintf('UnterminatedStrings', '--> n_units     : %g (all available units)\n', n_units);
    aux.cprintf('UnterminatedStrings', '--> duration_sec: %g ms\n', duration_sec);
    aux.cprintf('UnterminatedStrings', '    Reconstruction:\n');
    aux.cprintf('UnterminatedStrings', '--> causality   : %d\n', iscausal);
    aux.cprintf('UnterminatedStrings', '--> lags_ms     : %g ms\n', lags_ms);
    aux.cprintf('UnterminatedStrings', '--> binwidth    : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands     : %g\n', n_bands);
end



%% Correlation Coefficients between STIMULI
CCs = CC_stimuli( spec_st );


%%
switch data_type
    case 'SU'
        len_unit_list = height(tbl_strf);        
        unit_list = 1:len_unit_list;    % so sorting in this case

    case 'MUA'
        dummy = load(fullfile(load.path_to_data('Analysis'), 'idx_MUA_good_sorted_unit_thr(0-7).mat'));
        unit_list = dummy.unit_list;
        len_unit_list = length(unit_list);

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end



%%
un = 103;   % number of units to include

a = cellfun(@(X) mean(X), tbl_strf.CCest_dry(unit_list(1:un),:), 'UniformOutput', false);
b = cellfun(@(X) mean(X), tbl_strf.CCest_drr(unit_list(1:un),:), 'UniformOutput', false);

% CC of STRFs reconstructions
% CCest_dry = cell2mat( tbl_strf.CCest_dry(unit_list(1:un),:)  );
% CCest_drr = cell2mat( tbl_strf.CCest_drr(unit_list(1:un),:)  );
CCest_dry = cell2mat( a );
CCest_drr = cell2mat( b );


%%
figure( 0 + fignum );
% clf;
markersize = 32;
fontsize = 28;
fontsize_big = 42;
fontsize_bigger = 48;


Zest_dry = CCest_dry;
Zest_drr = CCest_drr;

mu_est_dry = median(Zest_dry);
mu_est_drr = median(Zest_drr);

% % Quantiles 
% qnt = [0.025, 0.975];
% q_est_dry = quantile(Zest_dry, qnt) - mu_est_dry;
% q_est_drr = quantile(Zest_drr, qnt) - mu_est_drr;
% fprintf('-- Using QUANTILEs!\n');

% SEM
q_est_dry = std(Zest_dry)/sqrt(size(Zest_dry,1));  q_est_dry(2,:) = q_est_dry;
q_est_drr = std(Zest_drr)/sqrt(size(Zest_dry,1));  q_est_drr(2,:) = q_est_drr;
fprintf('-- Using SEM\n');

switch upper(data_type)
    case 'SU'
        ax = subplot(1,2,1); 
        
    case 'MUA'
        ax = subplot(1,2,2);
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end

M = [mu_est_dry(:), mu_est_drr(:)];

C(1, 1, :) = aux.rpalette(sprintf('new%02d',1));
C(1, 2, :) = aux.rpalette(sprintf('new%02d',2));
h = superbar.superbar(1:n_drr, M, 'E', q_est_dry', 'BarFaceColor', C);
set(gca, 'XTick', 1:n_drr);
set(gca, 'XTickLabel', [drr.labels(drr.ordered)]);
xlim([0.5, 5.5]);

%{
C = [aux.rpalette(sprintf('new%02d',1))
    aux.rpalette(sprintf('new%02d',2))];
h = bar(M); 
set(gca, 'XTick', 1:n_drr, 'XTickLabel', drr.labels(drr.ordered));
h(1).FaceColor = aux.rpalette(sprintf('new%02d',1));
h(2).FaceColor = aux.rpalette(sprintf('new%02d',2));

dx = 0.2*h(1).BarWidth;
hold on
% ERROR BARs
er = errorbar((1:5)'-dx, mu_est_dry , q_est_dry(1,:), q_est_dry(2,:), 'CapSize',24);
er(2) = errorbar((1:5)'+dx, mu_est_drr , q_est_drr(1,:)', q_est_drr(2,:)', 'CapSize',24);
for k = 1:length(er)
    er(k).Color = [0 0 0];                            
    er(k).LineStyle = 'none'; 
    er(k).LineWidth = 4;
end
%}

% CCs
hold on
plth = plot(1:n_drr, CCs, 'ks:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k' );
hold off

set(gca, 'FontSize', fontsize);
ylim([0.0, 1.2]);
legend(h(1,:), {'$r_{dry}$ vs. $\hat{r}_{drr}$', '$r_{drr}$ vs. $\hat{r}_{drr}$'}, 'FontSize', fontsize_big);
xlabel('DRR', 'FontSize', fontsize_big);
title(sprintf('%d %s', size(Zest_dry,1), upper(data_type)), 'FontSize', fontsize_bigger);
set(gca, 'Box', 'On');

switch upper(data_type)
    case 'SU'
        set(gca, 'Position', [0.0970, 0.1100, 0.3928, 0.8038]);
        ylabel('$CC_{response}$', 'FontSize', fontsize_big);
        
    case 'MUA'
        set(gca, 'Position', [0.5279, 0.1100, 0.3928, 0.8038]);
        set(ax, 'YTickLabel', '');
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end




%% Two-ways ANOVA
figure( 10 + fignum );
clf;

y   = [Zest_dry(:); Zest_drr(:)]; 

% Fisher's z-transform
z = fisher_z_transform(y);

%
g1a = repmat(drr.labels(drr.ordered), size(Zest_dry,1), 1);
g1  = [g1a(:); g1a(:)]; 
g2a = repmat( {'$r_{dry}$ vs. $\hat{r}_{drr}$'}, size(Zest_dry));
g2b = repmat( {'$r_{drr}$ vs. $\hat{r}_{drr}$'}, size(Zest_drr));
g2  = [g2a(:); g2b(:)];
group = {g1, g2};
[p, tbl, stats] = anovan(z, group, ...
    'model','interaction',...
    'varnames',{'D/R', 'Unit Type'},...
    'display', false ...
);
%tbl = cell2table(tbl(2:end,:), 'VariableNames', tbl(1,:));
data.(data_type).anova2 = tbl;

fprintf('==========================================\n');
fprintf('- %s 2-ways ANOVA:\n', upper(data_type));
fprintf('==========================================\n');    
disp(tbl);

data.(data_type).multcompare = multcompare(stats, 'Dimension', [1,2]);
fprintf('- multcompare results:\n');
disp([(1:length(data.(data_type).multcompare))', data.(data_type).multcompare ]);

title(sprintf('%d %s', size(Zest_dry,1), upper(data_type)), 'FontSize', fontsize_bigger);
set(gca, 'FontSize', fontsize_big);

fprintf('\n');







