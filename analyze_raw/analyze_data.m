function st = analyze_data(Sft, H, binwidth, do_jackknife)
%
% function st = analyze_data(Sft, H, binwidth)
%
% 
% Description:
% Performs a series of analysis over the data. I run this function over all
% of my measurements for PSTH (spikes), MUA, or LFP data.
%

if ~iscell(Sft)
    Sft = {Sft};
end

%fs = 1/(1e-3*binwidth);     % (Hz) sampling rate if the spectrogram's (Sft) time axis 

% STRF
n_win   = 140;            % (smp) # STRF's lags
n_bands = size(Sft{1},1);  % # of frequency bands
n_splits= 20; % 12; 
n_smp   = size(Sft{1},2);  % # of frequency bands
assert(size(H,1) == n_smp);

% Analyze all available DRR cases
%drr = get_DRR_list_and_indices; 
n_drr     = size(H, 2);
train_drr = 1:n_drr;
test_drr  = 1:n_drr;    %train_drr;

% % Create the table that holds all the analyzed data
% variable_names = {'strf', 'r_est', 'mi', 'mtf'};
% n_var = length(variable_names);
% T     = cell2table(cell(1,n_var), 'VariableNames', variable_names);

n_lags    = (n_win - mod(n_win,2))/2 + 1;
st.n_splits = n_splits;
st.strf   = nan(n_bands, n_lags, n_drr);
st.r_test = nan(ceil(n_smp/n_splits), n_drr);
st.r_est  = nan(ceil(n_smp/n_splits), n_drr);
st.train_idx= nan(1, n_drr);
st.test_idx = nan(1, n_drr);

% Houtgast's MTF
%n_foct= 8;  % # of 1/3 octave band filters (TODO: pre-set; fix it!!!)
%st.mi  = nan(n_foct, 2, n_drr);
%st.mtf = nan(n_foct, n_drr);
st.JN       = cell(1, n_drr); 

%% Loop over all available DRR cases
valid_drr_cases = find(~isnan(sum(H,1)));

for k = valid_drr_cases
    [st.strf(:,:,k), st.r_est(:,k), ~, test_st, st.JN{k}] =...
        strfpkg.strf_train_and_test(Sft, H, binwidth, ...
        'n_win',        n_win, ...
        'n_splits',     n_splits, ...
        'train_drr',    train_drr(k),...
        'test_drr',     test_drr(k),...
        'do_jackknife', do_jackknife);
    
    if 1 == n_drr
        st.JN = st.JN{k};
    end
    
    st.r_test(:,k) = test_st.y;
    st.train_idx(k)= train_drr(k);
    st.test_idx(k) = test_drr(k);
    
    % Modulation transfer function
    %[st.mi(:,:,k), st.mtf(:,k), st.foct] = Houtgast(fs, Sft{k}, st.strf(:,:,k), H(:,k));

    % DEBUG
    %{
        figure(99);
        clf;
        plot(zca([st.r_test(:,k), st.r_est(:,k)]))
        cc = corrcoef([st.r_test(:,k), st.r_est(:,k)]);
        title(sprintf('CC: %g', cc(1,2)));
    %}
end







