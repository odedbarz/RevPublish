% plot_stats.m
%
% Plots error bars for Dry-vs-DDR reconstructions.


clc
fignum = 10;
verbose = 1;

setup_environment('../');


%% Plot properties
fontsize = 32;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(2*fontsize); %64;

markersize = 36;
linewidth = 5;



%% Analyze
data_type   = 'MUA';       % {'SU', MUA'}
switch data_type
    case 'SU'
        n_units = [10, 25, 50, 103];

    case 'MUA'
        n_units = [10, 25, 50, 103, 150, 241];    

    otherwise
        error('--> Unrecognized DATA_TYPE!');

end

drr = get_DRR_list_and_indices;

% Initialize the data
len_units = length(n_units);

dummy = nan(drr.n_drr, len_units);

tbl_mu = array2table( dummy );
tbl_mu.Properties.VariableNames = arrayfun(@(N) sprintf('CC_%d', N), n_units, 'uni', 0);
tbl_mu.Properties.RowNames = drr.labels(drr.ordered);

% Fast duplicate, MATLAB style
tbl_se = tbl_mu;
tbl_mu2 = tbl_mu;
tbl_se2 = tbl_mu;
tbl_mu3 = tbl_mu;
tbl_se3 = tbl_mu;


% Start loading...
fprintf('Loading the CC arrays...\n');

for nn = 1:len_units
    
    fn_path = '../_data/Analysis/';
    fn_name = sprintf('analyzed_cc_%s_units(%d).mat', data_type, n_units(nn));
    fn_fullfile = fullfile( fn_path, fn_name );
    
    fprintf('- Loading: %s\n', fn_name);

    % Load the data
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
    warning off
    data = load(fn_fullfile);
    warning on

    % Extract data: 
    CC      = data.CC;         % Sdry vs. Sest, averaged over time
    CC2     = data.CC2;        % Sdrr vs. Sest
    CC3     = data.CC3;        % Sdry vs. Sdrr

    CCt     = data.CCt;         % compares dry vs. est, as a function of time
    CCt2    = data.CCt2;        % compares drr vs. est
    CCt3    = data.CCt3;        % compares drr vs. est

    % Sdry-vs-Sest
    tbl_mu{:,nn} = squeeze( median(CC, 2) );
    tbl_se{:,nn} = squeeze( median(CC, 2) );
    
    % Sdrr-vs-Sest    
    tbl_mu2{:,nn} = squeeze( median(CC2, 2) );
    tbl_se2{:,nn} = squeeze( median(CC2, 2) );

    % Sdry-vs-Sdrr
    tbl_mu3{:,nn} = squeeze( median(CC3, 2) );
    tbl_se3{:,nn} = squeeze( median(CC3, 2) );
    
end
fprintf('Finished loading\n');







%% CC vs DRR (legend: units)
figure( 0 + fignum );
clf;

% CC
ax = gca;
plth = plot(tbl_mu.Variables, 's:', 'MarkerSize', markersize, 'linewidth', linewidth);
set(gca, 'FontSize', fontsize);

if strcmpi('MUA', data_type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end
ax.XTick = 1:height(tbl_mu);
ax.XTickLabel = tbl_mu.Row;
xlabel('DRR', 'FontSize', fontsize_big);
ylabel('CC', 'FontSize', fontsize_big);
aux.ctitle('Dry Stimulus vs. Reconstruction',...
    sprintf('(%d %s)', n_units(nn), data_type));

legend_str = arrayfun(@(N) sprintf('%d %s', N, data_type), n_units, 'uni', 0);
legend( legend_str, 'Location', 'southwest', 'FontSize', fontsize_bigger);





%% CC vs units (legend: DRR)
figure( 10 + fignum );
clf;

% CC
ax = gca;
plth = plot(tbl_mu.Variables', 's:', 'MarkerSize', markersize, 'linewidth', linewidth);
set(gca, 'FontSize', fontsize);

if strcmpi('MUA', data_type)
    arrayfun(@(I) set(plth(I), 'MarkerFaceColor', plth(I).Color), 1:length(plth));
end
ax.XTick = 1:size(tbl_mu, 2);

xlabeltick = arrayfun(@(N) sprintf('%d %s', N, data_type), n_units, 'uni', 0);
ax.XTickLabel = xlabeltick;
xlabel('Units', 'FontSize', fontsize_big);
ylabel('CC', 'FontSize', fontsize_big);
title('Dry Stimulus vs. Reconstruction', 'FontSize', fontsize_big);

ylim([0.5, 1.0]);
legend( tbl_mu.Row, 'Location', 'northeastoutside', 'FontSize', fontsize_bigger);



%%


