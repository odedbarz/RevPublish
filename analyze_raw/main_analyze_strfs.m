%
% main_ANALYZE_strfs.m
%
% Description:
% Run over the hole database and extracts and analyzes various parameters 
% of the STRFs.
%

clc
addpath('../');

verbose= 1;
fignum = [];
FigSetup;


%% DRR cases
drr = get_DRR_list_and_indices;
n_drr = drr.n_drr;


%% Load STRF data
data_type = 'MUA';     % {'PSTH', 'MUA'}

if ~exist(sprintf('%s_list', lower(data_type)), 'var')
    aux.cprintf('g', '-> Loading %s data...\n', data_type);
    
    % This file contains
    %     mua_list        1x437     cell                
    %     psth_list       1x437     cell                
    %     spec_list       1x2       cell                
    %     stim_list       1x2       cell                
    %     tbl_slc         1x20      table    
    fn.path = '../.data/';
    % fn.load = 'analysis_STRF_loopover_(14-May-2020).mat';
    switch upper(data_type)
        case 'PSTH'
            %fn.load = 'analysis_STRFs_PSTH_bw(2)ms_fbands(30)_win(10)ms_(15-May-2020)';
            fn.load = 'analysis_STRFs_PSTH_bw(10)ms_fbands(30)_win(25)ms_(16-May-2020).mat';
        case 'MUA'
            %fn.load = 'analysis_STRFs_MUA_bw(2)ms_fbands(30)_win(10)ms_(15-May-2020)';
            fn.load = 'analysis_STRFs_MUA_bw(10)ms_fbands(30)_win(25)ms_(16-May-2020).mat';
        otherwise
            error('-> Unidentified DATA_TYPE!');
    end
    load([fn.path, fn.load]);

else
    aux.cprintf('g', '-> Using current %s loaded data!\n', data_type);

end


%% Load responses
if ~exist('H', 'var')

    switch upper(data_type)
        case 'PSTH'
            %fn.load = 'data_PSTH_bw(2)_fbands(30)_win(20)ms_(07-May-2020).mat';
            fn.load = 'data_PSTH_bw(10)_fbands(30)_win(25)ms_(18-May-2020).mat';
        case 'MUA'
            %fn.load = 'data_MUA_bw(2)_fbands(30)_win(30)ms_(07-May-2020).mat';
            fn.load = 'data_MUA_bw(10)_fbands(30)_win(25)ms_(18-May-2020).mat';
        otherwise
            error('-> Unidentified DATA_TYPE!');
    end
    load([fn.path, fn.load]);

end



%% Set parameters & database
n_bands     = spec_list{1}.n_bands;     	% 
binwidth    = spec_list{1}.binwidth;     	% (ms)
win_size_ms = spec_list{1}.win_size_ms;    	% (ms)
f           = spec_list{1}.f;
assert( isequal(spec_list{1}.f, spec_list{2}.f), '-> All spectrograms has to have the same frequency bands!\n');


% Remove empty or NaN values
data_type_to_use = sprintf('%s_list', lower(data_type));
data        = eval( data_type_to_use );
nonempty_idx= cellfun(@(X) ~isempty(X), data);  % valid measurements' indices
% data        = data(nonempty_idx);
n_data      = length(data);

if verbose
    fprintf('-> data_type  : %s\n', data_type_to_use);
    fprintf('-> n_bands    : %d\n', n_bands);
    fprintf('-> binwidth   : %d (ms)\n', binwidth);
    fprintf('-> win_size_ms: %d (ms)\n', win_size_ms);
    fprintf('-> # of measurements: %d (%d are empty)\n', n_data, n_data-nnz(nonempty_idx));
end




%% Aggregate all data before showing it all
aux.cprintf('g', '-> Starting looping over the data...\n');

unit_rows  = nan(n_data, n_drr);        % keep the TBL_IMPALE row numbers of the units been used for the analysis
CF         = nan(n_data, n_drr);
cf_width   = nan(n_data, n_drr);
Tpeak      = nan(n_data, n_drr);
best_rate  = nan(n_data, n_drr);
best_scale = nan(n_data, n_drr);
power_strf = nan(n_data, 4); 
score      = nan(n_data, drr.n_drr);


for k = 1:n_data    
    if isempty(data{k})
        continue;
    end
    strf = data{k}.strf;

    pars = cell(1, size(strf,3));    
    for m = 1:n_drr       
        if any(any(isnan(strf(:,:,m))))
            continue;
        end
        pars{m} = strfpkg.analyze_strf(strf(:,:,m), f, binwidth, [] );
        
        unit_rows(k,m) = k;
        CF(k,m)        = pars{m}.cf.CF;                 % (kHz)
        cf_width(k,m)  = pars{m}.cf.width;              % (kHz)
        Tpeak(k,m)     = pars{m}.Tpeak;                 % (ms)
        best_rate(k,m) = pars{m}.mtf.args.best_rate;    % (Hz)
        best_scale(k,m)= pars{m}.mtf.args.best_scale;   % (Cycles/Octave)

    end
    
    % >>> Correlation between test & estimated response
    for q = 1:drr.n_drr
        dummy = corrcoef( data{k}.r_test(:,q), data{k}.r_est(:,q) );
        score(k, q) = dummy(1, 2);
    end
    
end


Qs = cf_width./CF;


%% Save the analysis results 
% %{
    'SAVE the STRF-PARS analysis!'    
    fn.save = sprintf('../.data/analyze_STRF-pars_%s_bw(%g)_fbands(%d)_win(%g)ms_(%s).mat',...
        data_type, binwidth, n_bands, win_size_ms, date);
    
    save(fn.save, '-v7.3', 'unit_rows', 'CF', 'cf_width', 'cf_width', 'Tpeak',...
        'best_rate', 'best_scale', 'spec_list');
%}


%% Sort the CFs and save them all into an excel file
[~, sorted_cfs_idx] = sort( CF(:,3) );




%% Analyze -- Distance between STRFs 
%{
dist_strf  = nan(n_data, 4);    % distance between DRY and other DRR conditions 
fdist      = @(X,Y) sum(abs(X(:)-Y(:)).^2); 
for k = 1:n_data    
    strf = data{k}.strf;
       
    % Get the DRY condition
    if ~isnan(sum(sum( strf(:,:,drr.dry) )))
        dry_idx = drr.dry;
    elseif isnan(sum(sum( strf(:,:,drr.dry) ))) && ~isnan(sum(sum( strf(:,:,drr.dry2) )))
        dry_idx = drr.dry2;
    else
        continue;
    end
    
    % >>> Compute distance between DRY & other DRR condiiotns
    for n = 1:4
        dist_strf(k,n) = fdist( strf(:,:,dry_idx), strf(:,:,drr.sortby(1+n)) );
    end
    
end


% Plot -- Change of STRFs as a function of DRR condition 
figure(100);
clf;
label_str = cell(1,4);
for n = 1:4
    label_str{n} = aux.ctitle('DRY vs.', sprintf('%s', drr.labels{drr.sortby(1+n)}) );
end
plot_dotbox(dist_strf, 'labels', label_str);
ylabel('$\Sigma |x - y|^2$');
title('Change of STRFs as a function of DRR condition');
%}


%% Analyze -- Compute STRFs power
%{
for k = 1:n_data    
    strf = data{k}.strf;

    for n = 1:drr.n_drr
        % Shechter_Depireux, 2007, Stability of spectro-temporal tuning over 
        % several seconds in primary auditory cortex of the awake ferret
        dummy_n = strf(:,:,n);
        power_strf(k,n)= sum( abs(dummy_n(:)).^2 );
    end
    

end
%}


%% Analyze -- Aggregate all test & estimates responses into big matrices
%{
% n_smp = size(data{1}.r_est, 1);
%n_smp = max(arrayfun(@(I) size(data{I}.r_est,1), 1:length(data)));  % takes time
n_smp = 300; % 334;    % for 40 ms duration
R     = nan(n_smp, n_data, drr.n_drr);
Rest  = nan(n_smp, n_data, drr.n_drr);
for k = 1:n_data    
    if size(data{k}.r_test,1) == n_smp
        R(:, k, :) = data{k}.r_test;
        Rest(:, k, :) = data{k}.r_est;
    end

end

%}




%% Plot results
% Define the measurement set to use

% I = 1:size(CF,1);
I = ~isnan(sum(CF(:,[1:5]),2));

nI = nnz(I);



%% Count invariant-CF neurons
analyze_strf_CF_invariant;


%%








%% Plot ** Compare Correlation Scores
figure(10);
clf;
n_hist = 25;
fontsize = 18;

ax = subplot(2,3,1);
[bincounts, bins] = histcounts(CF(:, drr.dry), n_hist);
bins = bins(1:end-1)+mean(diff(bins))/2;
bar(bins, bincounts);
xlabel('Frequency (kHz)');
title('CFs');
axis tight
ylabel('Count');
set(gca, 'FontSize', fontsize);

ax(2) = subplot(2,3,2);
[bincounts, bins] = histcounts(Tpeak(:, drr.dry), n_hist);
bins = bins(1:end-1)+mean(diff(bins))/2;
bar(bins, bincounts);
xlabel('Time (ms)');
title('Tpeak');
axis tight
set(gca, 'FontSize', fontsize);

ax(3) = subplot(2,3,3);
[bincounts, bins] = histcounts(Qs(:, drr.dry), n_hist);
bins = bins(1:end-1)+mean(diff(bins))/2;
bar(bins, bincounts);
xlabel('CF vs. $\Delta w$');
title('Q Factor');
axis tight
set(gca, 'FontSize', fontsize);

ax(4) = subplot(2,3,4);
[bincounts, bins] = histcounts(best_rate(:, drr.dry), n_hist);
bins = bins(1:end-1)+mean(diff(bins))/2;
bar(bins, bincounts);
xlabel('Rate (kHz)');
ylabel('Count');
title('Best Rate');
axis tight
set(gca, 'FontSize', fontsize);

ax(5) = subplot(2,3,5);
[bincounts, bins] = histcounts(best_scale(:, drr.dry), n_hist);
bins = bins(1:end-1)+mean(diff(bins))/2;
bar(bins, bincounts);
xlabel('Scale (Cycles/Octave)');
title('Best Scale');
axis tight
set(gca, 'FontSize', fontsize);


aux.abc(ax, [], 'northwestoutside');



%% Plot STRF vs MTF
% ne = [39, 403, 376];      % MUA
ne = [40, 96, 146];         % PSTH
figure(11);
clf

ax = zeros(1,length(ne));
f = 1e-3*spec_list{1}.f;     % (kHz)
for k = 1:length(ne)
    strf = data{ne(k)}.strf;
    ax(k) = subplot(3,1,k);
    h = strfpkg.plot_strf(strf(:,:,drr.dry), f, binwidth);    
    title('');
    if k == 1
        title('STRF');        
    elseif k<length(ne)
        set(ax(k), 'XTickLabel', '');
        xlabel(ax(k), '');
    end
    
    h = legend(h.Children(1), sprintf('CF: %.2f kHz', CF(ne(k),drr.dry)));
    h.Box = 'off';
end


%% Plot ** Compare CFs
figure(12);
clf;

% CFs
ax1 = zeros(1,4);
djitter = 0.25;
xjitter = djitter*rand(nI, 1);
yjitter = djitter*rand(nI, 1);
for k = 1:4
    ax1(k) = subplot(1,4,k);
    ii = drr.sortby(1+k);
    plot(CF(I,ii) + yjitter, CF(I,3) + xjitter, '.');
    xlabel(sprintf('$CF_{(%s)}$', drr.labels{ii}));
    %ylabel(sprintf('$CF_{(%s)}$', drr.labels{4}));
    %title('CF');
    axis tight square
    hold on
    plot(xlim, ylim, '--k')
    hold off
    dummy = corrcoef(CF(I,3), CF(I,ii));
    CC(k) = dummy(1,2);
    title(sprintf('$r^2: %.3f$', CC(k)));
end
axes(ax1(1));
ylabel(ax1(1), sprintf('$CF_{(%s)}$', drr.labels{3}) );




%% Compare Correlation Scores
figure(18);
clf;
ax1 = zeros(1,4);
djitter = 0.0;
xjitter = djitter*rand(n_data, 1);
yjitter = djitter*rand(n_data, 1);
CC = nan(1, 4);
for k = 1:4
    ax1(k) = subplot(1,4,k);
    ii = drr.sortby(1+k);
    plot(score(:,ii) + yjitter, score(:,3) + xjitter, '.');
    xlabel(sprintf('$CC_{(%s)}$', drr.labels{ii}));
    axis tight square
    hold on
    plot(xlim, ylim, '--k')
    hold off
    I = ~isnan(score(:,3)) & ~isnan(score(:,ii));    % remove nans
    dummy = corrcoef(score(I,3), score(I,ii));
    CC(k) = dummy(1,2);
    title(sprintf('$r^2: %.3f$', CC(k)));
end
axes(ax1(1));
ylabel(ax1(1), sprintf('$CC_{(%s)}$', drr.labels{3}) );



%%
figure(20);
clf;
f = 1e-3*spec_list{1}.f;
strf = data{39}.strf;

for n = 1:5
    ax = subplot(1,5,n);
    I = drr.sortby(n);
    strfpkg.plot_strf(strf(:,:,I), f, binwidth);
    title(sprintf('%s', drr.labels{ I }));

    if n > 1
        ylabel('');        
    end
    if n~=3
        xlabel('');        
    end

    set(gca, 'FontSize', fontsize);

end








