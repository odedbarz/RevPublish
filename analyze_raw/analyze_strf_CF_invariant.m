%
% analyze_strf_CF_invariant.m
%
% 

clc
addpath('../');

% Notes:
% * You need to run the **main_ANALYZE_strfs.m** analysis to have the
%   units' CFs,
%   OR, you can get load it

%   Name              Size            Bytes  Class     Attributes
% 
%   CF              437x6             20976  double              
%   Tpeak           437x6             20976  double              
%   best_rate       437x6             20976  double              
%   best_scale      437x6             20976  double              
%   cf_width        437x6             20976  double              
%   unit_rows       437x6             20976  double              
load('..\.data\analyze_STRF-pars_MUA_bw(10)_fbands(30)_win(50)ms_(01-Jun-2020).mat')


%%
drr = get_DRR_list_and_indices;


%% Claclulate the conditions for the invariance of the CFs
fignum = 50;

% Count invariant-CF neurons
ST = log2( 2.^([1, 3, 6, inf]/12) );   % (semitones) log2( 2^...) for clarity
n_dev = length(ST);
CF_ = CF./CF(:,drr.dry);
CF_ = CF_(:,drr.sortby(2:5));
not_nans = all(~isnan(CF_),2);
CF_ = CF_(not_nans,:);
legit_units = size(CF_,1);  % measurements that are all not nans
all_units = size(CF,1);
cnt = zeros(1, n_dev);
conds_cf = nan(legit_units, n_dev);

% Don't count the same neuron more than ones
used_neurons = false(legit_units,1);

for jj = 1:n_dev
    bnds = abs(log2(CF_)) <= ST(jj);
    conds_jj = all(bnds,2);
    if jj > 1
        % Exclude previous neurons
        conds_jj = conds_jj & ~used_neurons;
    end
    used_neurons = used_neurons | conds_jj;
    conds_cf(:,jj) = conds_jj;
    
    cnt(jj) = nnz( conds_jj );        
end
assert(legit_units == sum(cnt), '--> Something is wrong with the counting!');
assert(all( sum(conds_cf) == cnt ), '--> Something is wrong with the counting!');




%% Plot #1: 
%   plots histograms of CF count for deviation of CFs across reverberation
%   conditions
figure( fignum );
clf;
subplot(1,2,1);
bar(1:n_dev, cnt);
xticklabels = {'$0-1$', '$1-3$', '$3-6$', '$>6$'};
set(gca, 'XTickLabel', xticklabels);
xlabel('Semitones');
ylabel('Count');
title('CF-Invariance under Reverberation');
set(gca, 'FontSize', 20);
%
subplot(1,2,2);
hist(CF(not_nans,drr.dry), 25);
xlabel('Frequency (kHz)');
ylabel('Count');
title('Characteristic Frequencies');
set(gca, 'FontSize', 20);
legend(sprintf('%d Measurements', legit_units));




%% Plot #2
% select from the "invariant" neuron pull; column #1 are the "best"
% invariant units
neu        = find(conds_cf(:,1));  
neu        = neu([1, fix(end/2), end]);     % select 3 examples to show
unnorm_idx = find(not_nans);                % get the unnormalized CFs
n_neu      = length(neu);
fkHz       = 1e-3*spec_list{1}.f;

figure( fignum + 5);
clf;
ax = nan(n_neu,5);

for k = 1:n_neu
    I = unnorm_idx(neu(k));
    strf_ = data{I}.strf;

    for m = 1:5
        m_sorted = drr.sortby(m);
        ax(k,m) = subplot(n_neu,5,(k-1)*5 + m);
        strfpkg.plot_strf( strf_(:,:,m_sorted), fkHz, binwidth );
        
        % Set the CF text
        CF_km = CF(I, m_sorted);  % (kHz)
        text(0.4, 7, 1.2*max(strf_(:)), sprintf('CF: %.1f kHz', CF_km));
        
        if ~(n_neu == k) || ~(m == 1)
            set(gca, 'XTickLabel', '');
            set(gca, 'YTickLabel', '');
            xlabel('');
            ylabel('');
        end
        
        if 1 == k
            title(sprintf('%s', drr.labels{m_sorted}));
        else
            title('');
        end
    end   
end


ylabel(ax(2,1), 'Unit', 'FontSize', 24);



%% Check correlations between responses for CF-invariant units
% select from the "invariant" neuron pull; column #1 are the "best"
% invariant units
neu_cfinv  = find(conds_cf(:,1));  
unnorm_idx = find(not_nans);            % get the unnormalized CFs
n_neu      = length(neu_cfinv);
%fkHz       = 1e-3*spec_list{1}.f;
Mcc         = nan(n_neu, 4);    % xross-correlation between DRRs matrix; DRR x4: drr.sortby = [2,5,1,4] 
drrs        = [2,5,1,4];
    
for k = 1:n_neu
    Irow = unnorm_idx(neu_cfinv(k));
    I = find(tbl_slc.row == Irow);
    
    for m = 1:length(drrs)
        dummy = corrcoef(H(:,3,I), H(:,drrs(m),I));
        Mcc(k, m) = dummy(1,2);
    end
    
end


figure( fignum + 10 );
clf;
ax = nan(1, size(Mcc,2));
for k = 1:size(Mcc,2)
    ax(k) = subplot(1,size(Mcc,2),k);
    hist(Mcc(:,k), 10);
    xlim([0, 1]);
    ylim([0, 0.75*size(Mcc,1)]);
    title(drr.labels{drrs(k)});
    xlabel('CC');
end

ylabel(ax(1), 'Count');















