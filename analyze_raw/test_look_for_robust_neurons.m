%
% test_look_for_robust_neurons.m
%
%
% Description:
% Run over the hole database and search for neurons that are invariant to
% reverberation.
%

clc
addpath('../');

verbose= 1;
fignum = [];
FigSetup;


%% DRR cases
drr = get_DRR_list_and_indices;
n_drr = 5; % drr.n_drr;


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
            fn.load = 'data_PSTH_bw(2)_fbands(30)_win(20)ms_(07-May-2020).mat'; '** no STRF! **'
            %fn.load = 'analysis_STRFs_PSTH_bw(2)ms_fbands(30)_win(10)ms_(15-May-2020).mat';
            %fn.load = 'analysis_STRFs_PSTH_bw(10)ms_fbands(30)_win(25)ms_(16-May-2020).mat';
        case 'MUA'
            fn.load = 'data_MUA_bw(2)_fbands(30)_win(30)ms_(07-May-2020).mat'; '** no STRF! **'
            %fn.load = 'analysis_STRFs_MUA_bw(2)ms_fbands(30)_win(10)ms_(15-May-2020)';
            %fn.load = 'analysis_STRFs_MUA_bw(10)ms_fbands(30)_win(25)ms_(16-May-2020).mat';
        otherwise
            error('-> Unidentified DATA_TYPE!');
    end
    load([fn.path, fn.load]);

else
    aux.cprintf('g', '-> Using current %s loaded data!\n', data_type);

end




%% Set parameters & database
% n_bands     = spec_list{1}.n_bands;     	% 
% binwidth    = spec_list{1}.binwidth;     	% (ms)
% win_size_ms = spec_list{1}.win_size_ms;    	% (ms)
% f           = spec_list{1}.f;
% assert( isequal(spec_list{1}.f, spec_list{2}.f), '-> All spectrograms has to have the same frequency bands!\n');
% 
% 
% % Remove empty or NaN values
% data_type_to_use = sprintf('%s_list', lower(data_type));
% data        = eval( data_type_to_use );
% idx         = cellfun(@(X) ~isempty(X), data);  % valid measurements' indices
% data        = data(idx);
% n_data      = length(data);

n_bands     = spec_st.n_bands;     	% 
binwidth    = spec_st.binwidth;     	% (ms)
win_size_ms = spec_st.win_size_ms;   % (ms)
f           = spec_st.f;
n_data      = size(H,3);



if verbose
    %fprintf(verbose, '-> data_type  : %s\n', data_type_to_use);
    fprintf(verbose, '-> n_bands    : %d\n', n_bands);
    fprintf(verbose, '-> binwidth   : %d (ms)\n', binwidth);
    fprintf(verbose, '-> win_size_ms: %d (ms)\n', win_size_ms);
end

aux.cprintf('g', '-> # of measurements: %d\n', n_data);



%% Aggregate all data before showing it all
aux.cprintf('g', '-> Starting looping over the data...\n');

CF         = nan(n_data, n_drr);
% cf_width   = nan(n_data, n_drr);
% Tpeak      = nan(n_data, n_drr);
% best_rate  = nan(n_data, n_drr);
% best_scale = nan(n_data, n_drr);
% dist_strf  = nan(n_data, 4);    % distance between DRY and other DRR conditions 
% fdist      = @(X,Y) norm(X(:)-Y(:),2); 
% power_strf = nan(n_data, 4); 
% score.est = nan(n_data, drr.n_drr);
% score.err = nan(n_data, drr.n_drr);

avg_rate1  = nan(n_data, n_drr);
avg_rate2 = nan(n_data, n_drr);
std_rate1  = nan(n_data, n_drr);
std_rate2  = nan(n_data, n_drr);

n_phoneme = 400/binwidth;
i_phoneme1 = 1790 + (1:n_phoneme);
i_phoneme2 = 15370 + (1:n_phoneme);

for k = 1:n_data   
    %{ 
    strf = data{k}.strf;
    
    pars = cell(1, size(strf,3));    
    for m = 1:n_drr       
        if any(any(isnan(strf(:,:,m))))
            continue;
        end
        pars{m} = strfpkg.analyze_strf(strf(:,:,m), f, binwidth, fignum);
        CF(k,m) = pars{m}.cf.CF;     % (kHz)
        
        % >>> Correlation between test & estimated response
        dummy = corrcoef( data{k}.r_test(:,m), data{k}.r_est(:,m) );
        score.est(k, m) = dummy(1, 2);
        
        am = strf(:,:,drr.dry);
        am = (am - mean(am(:)))/std(am(:));
        bm = strf(:,:,m);
        bm = (bm - mean(bm(:)))/std(bm(:));        
        d_strf = am(:) - bm(:);
        score.err(k, m) = rms( d_strf );
        
        avg_rate(k,m) = mean(data{k}.r_test( i_phoneme ,m));
        std_rate(k,m) = std(data{k}.r_test( i_phoneme ,m));
    end
    %}
    
    avg_rate1(k,:) = mean( H( i_phoneme1, 1:n_drr, k));
    avg_rate2(k,:) = mean( H( i_phoneme2, 1:n_drr, k));
  
    std_rate1(k,:) = std( H( i_phoneme1, 1:n_drr, k));
    std_rate2(k,:) = std( H( i_phoneme2, 1:n_drr, k));
    
end


%%

figure(11);
clf;
for ii = 1:n_data
    %ax = subplot(1,2,1);
    hold on    
%     plot(ii*ones(1,n_drr), avg_rate1(ii,:)', '.');
    plot(ii*ones(1,n_drr), std_rate1(ii,:)', '.');
    %plot(ii*ones(1,n_drr), std_rate1(ii,:)'./avg_rate1(ii,:)', '.');
    
    %ax(2) = subplot(1,2,2);  
    %hold on
%     plot(ii*ones(1,n_drr), avg_rate2(ii,:)', 's', 'MarkerSize', 9);
    plot(ii*ones(1,n_drr), std_rate2(ii,:)', 's', 'MarkerSize', 9);
    %plot(ii*ones(1,n_drr), std_rate2(ii,:)'./avg_rate2(ii,:)', 's', 'MarkerSize', 9);
end

hold off
% linkaxes(ax);












