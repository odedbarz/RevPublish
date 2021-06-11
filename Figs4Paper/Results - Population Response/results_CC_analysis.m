%
% results_CC_analysis.m
%

clc
fignum = 10;
verbose = 1;

addpath('../../');
setup_environment;


%% Load the data
%              CC: [7200×5 double] 	% compares dry vs. est
%             CC2: [7200×5 double]	% compares drr vs. est
%             CCt: [7200×5 double]  % compares dry vs. est along the time domain
%            CCt2: [7200×5 double]  % compares drr vs. est along the time domain
%             PPt: [7200×5 double]
%            PPt2: [7200×5 double]
%     patch_width: 1
%         spec_st: [1×1 struct]
%          splits: [1×1 struct]
%         stim_st: [1×1 struct]
%        tbl_data: [103×20 table]
data_type_list = {'SU', 'MUA'};
FN_TEMPLATE = 'analyzed_%s_CC_(07-Jun-2021)_units(%d).mat';
FN_PATH = load.path_to_data('Analysis');

% Available units:
% 'SU' : [10, 25, 50, 103] 
% 'MUA': [10, 25, 50, 103, 150, 241]   
% units = [10, 25, 50, 103];    % units to load
units = [10, 25, 50, 100]; %, 150, 241]; 
drr = get_DRR_list_and_indices;
n_drr = drr.n_drr;

% T.units = units;
% T.su.mu = array2table(nan(n_drr, length(units)));
% T.su.mu.Properties.RowNames = drr.labels(drr.ordered);
% T.su.mu.Properties.VariableNames = arrayfun(@(N) sprintf('unit%d', N), units, 'UniformOutput', false);
% T2add.su.mu = array2table(nan(n_drr, length(units)));
% 
% T.su.SE = T.su.mu;  % duplicate for SE
% T.mua = T.su;       % duplicate for MUA

T = array2table(nan(n_drr, length(units)));
T.Properties.VariableNames = arrayfun(@(N) sprintf('unit%d', N), units, 'UniformOutput', false);

tbl_spk = array2table( cell(n_drr, length(units)) );
tbl_spk.Properties.VariableNames = arrayfun(@(N) sprintf('unit%d_spk', N), units, 'UniformOutput', false);

tbl_drr_labels = table(drr.labels(drr.ordered)', 'VariableNames', {'DRR'});
tbl_type = table(cell(n_drr,1), 'VariableNames', {'type'});

Tmu_0 = [tbl_drr_labels, tbl_type, T, tbl_spk];
Tse_0 = [tbl_drr_labels, tbl_type, T];
clear T     % just a template for the other tables

% *** Sdry vs. Sest ***
Tmu = [];
Tse = [];

% *** Sdrr vs. Sest ***
Tmu2 = [];
Tse2 = [];

n_speakers = 12;


for data_type = data_type_list
    data_type = data_type{:};
    Tmu_k = Tmu_0;
    Tse_k = Tse_0;    
    Tmu2_k = Tmu_0;
    Tse2_k = Tse_0;        
    for un = units
        if strcmpi('SU', data_type) && un > 103
            col_name = sprintf('unit%d', un);
            Tmu_k.([col_name, '_spk']) = nan(n_drr, n_speakers);
            Tmu2_k.([col_name, '_spk']) = nan(n_drr, n_speakers);
            continue;
        end
        
        % Load the data:
        fn_name = sprintf(FN_TEMPLATE, data_type, un);
        fn_fullfile = fullfile( FN_PATH, fn_name );
        warning off
        data = load(fn_fullfile);
        warning on

        aux.cprintf('String', '\nLoading analyzed data...\n');
        aux.cprintf('String', 'data_type: %s\n', upper(data_type));
        aux.cprintf('String', 'n_units  : %d\n', un);
        aux.cprintf('String', 'filename : %s\n', fn_name);
        aux.cprintf('String', 'path     : %s\n', FN_PATH);
        
        % *** Sdry vs. Sest ***
        [Tmu_k, Tse_k] = S_vs_Sest(data.tbl.CC.Variables, data_type, un, Tmu_k, Tse_k);
        
        % *** Sdrr vs. Sdrr ***
        [Tmu2_k, Tse2_k] = S_vs_Sest(data.tbl.CC2.Variables, data_type, un, Tmu2_k, Tse2_k);
        
    end
    
    Tmu = [Tmu; Tmu_k];    
    Tse = [Tse; Tse_k];
    
    Tmu2 = [Tmu2; Tmu2_k];    
    Tse2 = [Tse2; Tse2_k];
    
end



%% Summary data for plotting
idx.su = contains(Tmu.type, 'SU');
idx.mua = contains(Tmu.type, 'MUA');

% Set the number of units to use for the comparison
'================================='
num_units_to_analyze = 100
'================================='

x_bias = 0.25;  

cc_units = sprintf('unit%d', num_units_to_analyze);


clear su
su.CCmu_units   = cc_units;
su.CCmu_info    = 'DRRs x SPEAKERS';
su.CCunits_info = sprintf('DRRs x UNITS (%d-%d)', units(1), units(end));
su.CCunits      = nan(n_drr, length(units));

clear mua
mua.CCmu_units  = cc_units;
mua.CCmu_info   = 'DRRs x SPEAKERS';
mua.CCunits_info= 'DRRs x UNITS (10-241)';
mua.CCunits     = nan(n_drr, length(units));

% Sdry vs Sest
su.CCspk  = Tmu(idx.('su'), [cc_units, '_spk']).Variables;   % table to array
mua.CCspk = Tmu(idx.('mua'), [cc_units, '_spk']).Variables;   % table to array

% Sdrr vs Sest
su.CCspk2  = Tmu2(idx.('su'), [cc_units, '_spk']).Variables;   % table to array
mua.CCspk2 = Tmu2(idx.('mua'), [cc_units, '_spk']).Variables;   % table to array


% su.CCunits_info = sprintf('DRRs x UNITS (%d-%d)', units(1), units(end));
% su.CCunits = nan(n_drr, length(units));
% 
% mua.CCunits_info = 'DRRs x UNITS (10-241)';
% mua.CCunits = nan(n_drr, length(units));

for k = 1:length(units)
    col_name = sprintf('unit%d', units(k));
    
    % Sdry vs Sest
    su.CCunits(:,k) = median( Tmu(idx.('su'), [col_name, '_spk']).Variables, 2 );
    mua.CCunits(:,k) = median( Tmu(idx.('mua'), [col_name, '_spk']).Variables, 2 );
    
    % Sdrr vs Sest    
    su.CCunits2(:,k) = median( Tmu2(idx.('su'), [col_name, '_spk']).Variables, 2 );
    mua.CCunits2(:,k) = median( Tmu2(idx.('mua'), [col_name, '_spk']).Variables, 2 );
    
end




%% Correlation Coefficients between STIMULI
CCs = CC_stimuli( data.spec_st );




%% SU-vs-DRR & MUA-vs-DRR side by side
figure(0+fignum);
clf;

linewidth = 5;
markersize = 25;
fontsize = 28;
fontsize_big = 38;
fontsize_bigger = 45;


% drr_idx = drr.sortby(1:n_drr);

% *** SU ***
ax = subplot(1,2,1);
plth = plot(su.CCunits(:,1:4), 's:', 'MarkerSize', markersize, 'LineWidth', linewidth);

% Stimulus-only performance (CCs)
hold on
plot(CCs, 'sk:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k', 'DisplayName', 'CC (Stim Dry-to-DRR)');
hold off

set(gca, 'FontSize', fontsize);
set(gca, 'XTick', 1:n_drr);
set(gca, 'XTickLabel', drr.labels(drr.sortby));
xlabel(ax(1), 'DRR', 'FontSize', fontsize_big );
ylabel('$CC (S_{dry}-to-\hat{S})$', 'FontSize', fontsize_big );
title('SU Reconstructions', 'FontSize', fontsize_bigger );
units_list = arrayfun(@(n) sprintf('%d SUs', n), units(1:4), 'UniformOutput', 0);
legend_list = [units_list, 'CC (Stim Dry-to-DRR)'];
legend_h = legend(legend_list, 'Location', 'southwest');
legend_h.FontSize = fix(0.8*fontsize);


% *** MUA ***
ax(2) = subplot(1,2,2);
plth = plot(mua.CCunits, 's:', 'MarkerSize', markersize, 'LineWidth', linewidth);
colors = get(plth,'Color');
if 1 < size(colors,1)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', colors{I}), 1:length(colors) );
else
    set(plth, 'MarkerFaceColor', colors);
end

% Stimulus-only performance (CCs)
hold on
plot(CCs, 'sk:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k');
hold off

set(gca, 'FontSize', fontsize);
set(gca, 'XTick', 1:n_drr);
set(gca, 'XTickLabel', drr.labels(drr.sortby));
xlabel(ax(2), 'DRR', 'FontSize', fontsize_big );
title('MUA Reconstructions', 'FontSize', fontsize_bigger );
units_list = arrayfun(@(n) sprintf('%d SUs', n), units, 'UniformOutput', 0);
legend_list = [units_list, 'CC (Stim Dry-to-DRR)'];
legend_h = legend(legend_list, 'Location', 'southwest');
legend_h.FontSize = fix(0.8*fontsize);

linkaxes(ax);
ylim([0.3, 1.0]);
xlim(ax, [0.7, 5.3]);
set(ax(2), 'YTickLabel', '');
ylabel(ax(2), '');

ax(1).Position = [0.1300    0.1100    0.3823    0.8038];
ax(2).Position = [0.5225    0.1100    0.3825    0.8038];





%% BOXPLOT
figure(10+fignum);
clf;

% su.grp  = -1 + 2*(1:n_drr) + x_bias*(-1);
% mua.grp = -1 + 2*(1:n_drr) + x_bias;
xspace = 4;
su.grp  = 0:xspace:5*xspace-1;
mua.grp = 1:xspace:5*xspace;

su_hnd = boxplot(su.CCspk', 'positions', su.grp, 'labels', su.grp, 'Colors', 'r');
hold on
mua_hnd = boxplot(mua.CCspk', 'positions', mua.grp, 'labels', mua.grp);
hold off
ylim([0.4, 1.0]);
set(gca, 'XTickLabel', drr.labels(drr.ordered), 'FontSize', fontsize);
assert( isequal(su.CCmu_units, mua.CCmu_units) );
title_str = sprintf('SU vs. MUA (%d units)', sscanf(mua.CCmu_units, 'unit%d'));
title( title_str, 'FontSize', fontsize_big);
legend([su_hnd(end-2), mua_hnd(end-2)], {'SU', 'MUA'}, 'FontSize', fontsize);
ylabel('CC', 'FontSize', fontsize_big);
xlabel('DRR', 'FontSize', fontsize_big);



%% Wilcoxon signed rank test between SU & MUA 
hnull = nan(1, n_drr);
pv = nan(1, n_drr);
for k = 1:n_drr
    [pv(k), hnull(k)] = signrank(su.CCspk(k,:), mua.CCspk(k,:));
end

fprintf('\n Wilcoxon signed rank test between SU & MUA\n')
array2table([pv; hnull], 'VariableNames', {'dry', 'dB9_4', 'dB4_8', 'dBm2_5', 'dBm8_2'} )




%% Plot BARS of SU & MUA 
figure(30 + fignum);
clf;

% Standard error
SE_fun = @(x) std(x,[],2)/sqrt(size(x,2));

% *** Sdry vs. Sest
ax = subplot(1,2,1);
su.median = [median(su.CCspk,2), median(su.CCspk2,2)];
su.SE     = [SE_fun(su.CCspk), SE_fun(su.CCspk2)];
title_str = sprintf('SU (%d units)', sscanf(mua.CCmu_units, 'unit%d'));
[h, plth, er] = plot_CC_S_vs_Sest(su.median, su.SE, CCs, title_str);
legend(ax(1), 'Off');

% *** Sdry vs. Sest
ax(2) = subplot(1,2,2);
mua.median = [median(mua.CCspk,2), median(mua.CCspk2,2)];
mua.SE     = [SE_fun(mua.CCspk), SE_fun(mua.CCspk2)];
title_str  = sprintf('MUA (%d units)', sscanf(mua.CCmu_units, 'unit%d'));
[h, plth, er] = plot_CC_S_vs_Sest(mua.median, mua.SE, CCs, title_str);

set(ax(2), 'YTickLabel', '');
ylabel(ax(2), '');
ax(1).Position = [0.1300    0.1100    0.3823    0.8038];
ax(2).Position = [0.5225    0.1100    0.3825    0.8038];

linkaxes(ax);




%% Wilcoxon signed rank test between SU & MUA 
tbl_pv.su = array2table(nan(2, n_drr), 'VariableNames', {'dry', 'dB9_4', 'dB4_8', 'dBm2_5', 'dBm8_2'} );
tbl_pv.mua = array2table(nan(2, n_drr), 'VariableNames', {'dry', 'dB9_4', 'dB4_8', 'dBm2_5', 'dBm8_2'} );

for k = 1:n_drr
    [tbl_pv.su{1,k}, tbl_pv.su{2,k}] = signrank(su.CCspk(k,:), su.CCspk2(k,:));
    [tbl_pv.mua{1,k}, tbl_pv.mua{2,k}] = signrank(su.CCspk(k,:), su.CCspk2(k,:));
end

fprintf('\n Wilcoxon signed rank test between SU (reverberant vs. dry speech)\n')
tbl_pv.su

fprintf('\n Wilcoxon signed rank test between MUA (reverberant vs. dry speech)\n')
tbl_pv.mua



%% Plot BARS of SU & MUA with ** SUPERPLOT **
figure(40 + fignum);
clf;

import superbar.*

% Standard error
SE_fun = @(x) std(x,[],2)/sqrt(size(x,2));

C = [aux.rpalette(sprintf('new%02d',1))
    aux.rpalette(sprintf('new%02d',2))];
xjupms = 3;
xbias = repmat([1 2], n_drr, 1) + (0:n_drr-1)'*xjupms;
su.median = [median(su.CCspk,2), median(su.CCspk2,2)];
su.SE     = [SE_fun(su.CCspk), SE_fun(su.CCspk2)];

% *** Sdry vs. Sest
ax = subplot(1,2,1);
for k = 1:n_drr
    if 0 == tbl_pv.su{2,k}
        superbar(xbias(k,:), su.median(k,:), 'E', su.SE(k,:), 'BarFaceColor', C);
    else
        P = nan(2);
        P(1,2) = tbl_pv.su{1,k};
        P(2,1) = P(1,2);
        superbar(xbias(k,:), su.median(k,:), 'E', su.SE(k,:), 'BarFaceColor', C,...
            'P', P);
    end
    hold on
end
plth = plot( mean(xbias,2), CCs, 'ks:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k' );
hold off
set(gca, 'FontSize', fontsize);
title_str = sprintf('SU (%d units)', sscanf(mua.CCmu_units, 'unit%d'));
title(title_str, 'FontSize', fontsize_bigger);
ylabel('CC', 'FontSize', fontsize_big);
xlabel('DRR', 'FontSize', fontsize_big);
set(gca, 'XTick', mean(xbias,2), 'XTickLabel', drr.labels(drr.ordered));
set(gca, 'Box', 'On');

% *** Sdry vs. Sest
ax(2) = subplot(1,2,2);
for k = 1:n_drr
    if 0 == tbl_pv.su{2,k}
        h_mua = superbar(xbias(k,:), mua.median(k,:), 'E', mua.SE(k,:), 'BarFaceColor', C);
    else
        P = nan(2);
        P(1,2) = tbl_pv.mua{1,k};
        P(2,1) = P(1,2);
        superbar(xbias(k,:), mua.median(k,:), 'E', mua.SE(k,:), 'BarFaceColor', C,...
            'P', P);
    end
    hold on
end
plth = plot( mean(xbias,2), CCs, 'ks:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k' );
hold off
set(gca, 'FontSize', fontsize);
title_str  = sprintf('MUA (%d units)', sscanf(mua.CCmu_units, 'unit%d'));
title(title_str, 'FontSize', fontsize_bigger);
xlabel('DRR', 'FontSize', fontsize_big);
set(gca, 'XTick', mean(xbias,2), 'XTickLabel', drr.labels(drr.ordered));
legend(h_mua, {'$S_{dry}$ vs. $\hat{S}_{drr}$', '$S_{drr}$ vs. $\hat{S}_{drr}$'},...
    'FontSize', fontsize_bigger);
set(gca, 'Box', 'On');
set(ax(2), 'YTickLabel', '');
ylabel(ax(2), '');
ax(1).Position = [0.1300    0.1100    0.3823    0.8038];
ax(2).Position = [0.5225    0.1100    0.3825    0.8038];

linkaxes(ax);
ylim([0, 1.3]);





