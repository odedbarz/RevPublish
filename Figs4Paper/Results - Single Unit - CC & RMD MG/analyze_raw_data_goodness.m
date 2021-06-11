%
% analyze_raw_data_goodness.m
%
% Description:
% Analyze the raw data file,
%   stats_raw(MUA)_(19-Feb-2021)_BW(5)ms_duration(36)sec_units(241).mat
%
% To create this file, analyze_raw/run main_raw_MUA_stats.m
%


clc
fignum = 10;
verbose = 1;

setup_environment('../../');

drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);
drr_ordered_labels = drr.labels(drr.ordered);



%% Load RAW data
%
%   Name           Size             Bytes  Class     Attributes
%   stats          1x1             270282  struct              
%   tbl_MUA      241x21            383652  table               
fn.path = load.path_to_data('Stats');
fn.file = 'stats_raw(MUA)_(02-Jun-2021)_BW(5)ms_duration(36)sec_units(241).mat';
load( fullfile(fn.path, fn.file) );



%%
figure(10+fignum);
clf;

fontsize = 28;	
fontsize_big = fix(1.5*fontsize); 
%fontsize_bigger = fix(3*fontsize); 
markersize = 28;

idx_drr = [1, 2, 5];
threshold = 0.7;

[sorted_units, idx_good_sorted_unit] = sort(stats.CC.median(idx_drr(1),:), 'descend');
idx_thr_unit   = find(stats.CC.median(idx_drr(1), idx_good_sorted_unit) >= threshold, 1, 'last');
good_units      = @(n) stats.CC.median(idx_drr(n), idx_good_sorted_unit(1:idx_thr_unit))';
not_good_units  = @(n) stats.CC.median(idx_drr(n), idx_good_sorted_unit(idx_thr_unit+1:end))';

spks_and_median = tbl_MUA.SPK .* stats.CC.median(idx_drr(1),:)';
idx_good_and_spks = find( spks_and_median >= threshold );
good_units_and_spk= @(n) stats.CC.median(idx_drr(n), idx_good_and_spks)';

plth = plot(good_units(1), [good_units(2), good_units(3)], '.');
set(plth, 'MarkerSize', fix(1.5*markersize));
axis square, xlim([0,1]), ylim([0, 1]);
hold on
plth2 = plot(not_good_units(1), [not_good_units(2), not_good_units(3)], 's');
set(plth2(1), 'MarkerSize', fix(0.45*markersize), 'Marker', 's', 'Color', aux.rpalette(1));
set(plth2(2), 'MarkerSize', fix(0.45*markersize), 'Marker', 's', 'Color', aux.rpalette(2));

plth2 = plot(good_units_and_spk(1), [good_units_and_spk(2), good_units_and_spk(3)], 'ok');

plot([0,1], [0,1], '--k');
plot(threshold*[1 1], [0,1], ':k');
hold off

legend_str1 = sprintf('median(%s)', drr_ordered_labels{idx_drr(2)});
legend_str2 = sprintf('median(%s)', drr_ordered_labels{idx_drr(3)});
legend({legend_str1, legend_str2}, 'Location', 'northwest');
xlabel('Within Sample Median (Dry)', 'fontsize', fontsize_big);
ylabel('Within Sample Median (DRR)', 'fontsize', fontsize_big);
set(gca, 'FontSize', fontsize);
title('Median Within Raw Samples (MUA)', 'fontsize', fontsize_big);

fprintf('\n');
fprintf('- Number of above-threshold units: %g\n', idx_thr_unit);
fprintf('- Number of Spike measurements AND above-threshold: %g\n', length(idx_good_and_spks));



%%
figure(15+fignum);
clf;

fontsize = 28;	
fontsize_big = fix(1.5*fontsize); 
fontsize_bigger = fix(3*fontsize); 

idx_drr = [1, 2, 5];

good_units_std     = @(n) stats.CC.std(idx_drr(n), idx_good_sorted_unit(1:idx_thr_unit))';
not_good_units_std = @(n) stats.CC.std(idx_drr(n), idx_good_sorted_unit(idx_thr_unit+1:end))';

% plth = plot(stats.std(1,:), [stats.std(idx_drr(1),:)', stats.std(idx_drr(2),:)'], '.');
plth = plot(good_units_std(1), [good_units_std(2), good_units_std(3)], '.');
set(plth, 'MarkerSize', markersize)
hold on
plth2 = plot(not_good_units_std(1), [not_good_units_std(2), not_good_units_std(3)], 's');
set(plth2(1), 'MarkerSize', fix(0.45*markersize), 'Marker', 's', 'Color', aux.rpalette(1));
set(plth2(2), 'MarkerSize', fix(0.45*markersize), 'Marker', 's', 'Color', aux.rpalette(2));
plot([0,0.3], [0,0.3], '--k');
hold off

legend_str1 = sprintf('median(%s)', drr_ordered_labels{idx_drr(2)});
legend_str2 = sprintf('median(%s)', drr_ordered_labels{idx_drr(3)});
legend({legend_str1, legend_str2}, 'Location', 'northwest');
xlabel('Within Sample STD (Dry)', 'fontsize', fontsize_big);
ylabel('Within Sample STD (DRR)', 'fontsize', fontsize_big);
set(gca, 'FontSize', fontsize);
title('STD Within Raw Samples (MUA)', 'fontsize', fontsize_big);

axis equal







