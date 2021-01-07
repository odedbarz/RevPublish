function psth_st = PSTH(S, duration_ms, binwidth, spikechan, syncchan, use_hcompute, fignum)
%
%   psth_st = psth(t, ch, duration_ms, binwidth, spikechan, syncchan, use_hcompute, fignum)
%
% Description:
% Transforms spike trials to PSTH(s).
%
% h holds the spike counts at each bin.
%
% if USE_HCOMPUTE == true compute also all the PSTHs for all trials, before
% averaging.

t = S.t;
ch= S.ch; 

if 3 > nargin || isempty(binwidth)
    binwidth  = 1;  % (ms)
end

if 5 > nargin || isempty(syncchan)
    syncchan  = 0;  % (ms)
end

if 4 > nargin || isempty(spikechan)
    %spikechan  = 1;  % (ms)    
    ch_list = unique( [ch{1}] );
    ch_list = ch_list(ch_list ~= syncchan);
    assert(~isempty(ch_list), '--> ERROR at [PSTH.m]: there are no spikes in these measurements!!');
    assert(length(ch_list)==1, '--> ERROR at [PSTH.m]: there are TOO MANY spikes in these measurements -- you need to CHOOSE ONE!!');
end

if 6 > nargin || isempty(use_hcompute)
    use_hcompute = 1;
end

if 7 > nargin || isempty(fignum)
    fignum = [];
end


%%
histSize = ceil(duration_ms/binwidth);

if use_hcompute
    h = set( histogram, ...
        'Type', 'PST', ...          % histogram's type
        'Size', [1, histSize], ...
        'BinWidths', binwidth, ...
        'Offsets', 0, ...
        'SyncChan', syncchan, ...
        'SpikeChan', spikechan ...
    );

    hc         = cellfun(@(T,CH) hcompute(h,T,CH), t, ch, 'UniformOutput', false);      % histograms's handles
    h          = cellfun(@(HC) double(HC)', hc, 'UniformOutput', false);              	% the PSTs cell(s)   
    bins       = get(hc{1}, 'BinTimes');
    nsync      = cellfun(@(HC) get(HC,'SyncCount'), hc, 'UniformOutput', true);            
    spikecount = cellfun(@(HC) get(HC,'SpikeCount'), hc, 'UniformOutput', true);        

else
    bins       = (0 : histSize) * binwidth;
    tps        = cellfun(@(T,CH) spiketools.pstimes(T,CH), t, ch, 'UniformOutput', false);
    kspike     = cellfun(@(CH) find(CH == spikechan), ch, 'UniformOutput', false);
    hc         = cellfun(@(TPS, KS) histc(TPS(KS), bins), tps, kspike, 'UniformOutput', false);
    h          = cellfun(@(HC) HC(1:end-1), hc, 'UniformOutput', false);    % discard last, overflow bin
    nsync      = cellfun(@(CH) nnz(0==CH), ch, 'UniformOutput', true);            
    spikecount = cellfun(@(KS) nnz(KS), kspike, 'UniformOutput', true);   
    
    % Remove the extra value
    bins = bins(1:end-1);

    % Also compute PSTHs for all trials separately
    trials      = cellfun(@(CH) cumsum(CH == syncchan), ch, 'UniformOutput', false);
    n_psths     = length(trials(:));
    H_trials    = cell(size(ch));
    nspk_trials= cell(size(ch));
    for p = 1:n_psths
        n_trials = max(trials{p});
        H_trials{p} = nan(size(h{p},1), n_trials);
        nspk_trials{p} = nan(1, n_trials);    % # of spikes in each trial
        for q = 1:n_trials
            tps_idx_pq = (ch{p} == spikechan) & (q == trials{p});
            tps_k = tps{p}(tps_idx_pq);
            H_trials{p}(:,q) = histc(tps_k, bins);
            nspk_trials{p}(q) = numel(tps_k);
        end
    end
    
    % Remove empty PSTH entries
    %H_trials = H_trials(nsync(:) > 0);

end


%% Remove empty PSTH entries
%h = h(nsync(:) > 0);


%% Normalize to units of (spikes/sec)
psth = cell(size(h));
% binwidth_sec = units.ms2sec(binwidth);
for k = 1:length(psth(:))
    if 0 < nsync(k)
        % Normalize to units of (spikes/sec)
        % psth{kk} = h{kk}/(binwidth_sec*n_sync(kk));    
        
        % Normalize by # of measurements
        psth{k} = h{k}/nsync(k);    
    end
end


%% REVERBERATION project -- assign the measurements in the appropriate order!
available_meas_idx = Impale_entries_to_indices(S);
n_drr     = length(available_meas_idx);
n_smp     = size(psth{1}, 1);
psth_st.H = nan(n_smp, n_drr);

% Exception: usually, I performed a measurement without the 3.0m\100% entry
% to save time (it's exactly 1.5m\100%)
if 5 == length(psth)
    available_meas_idx(6) = false;
end

psth_st.H(:,available_meas_idx) = [psth{:}];

% for k = 1:nnz(available_meas_idx)
%     if available_meas_idx(k)
%         psth_st.H(:, ) = psth{k};    % convert into a matrix
% 
% end



%% Set the output
psth_st.type         = 'psth';
psth_st.duration_ms  = duration_ms;
psth_st.binwidth     = binwidth;
psth_st.fs           = 1/units.ms2sec(psth_st.binwidth);          % (Hz)    
psth_st.spikechan    = spikechan;
psth_st.syncchan     = syncchan; 
psth_st.use_hcompute = use_hcompute;
psth_st.bins         = bins;
psth_st.nsync        = nsync;
psth_st.spikecount   = spikecount;

psth_st.available_meas_idx = available_meas_idx;
% pst_st.psth         = psth;
% psth_st.H            = [psth{:}];    % convert into a matrix

if ~use_hcompute
    psth_st.H_trials     = H_trials;
    psth_st.nspk_trials  = nspk_trials;
end
    
%% Plot  
if isempty(fignum), return; end
figure(fignum);
clf;
p1 = double(psth{1});
plot(bins, p1)

if exist('H_trials', 'var')
    plot(bins, h{1}, '-', bins, H_trials{1}(:,1), '-');
    legend('1 trial', 'PSTH');
    xlabel('Time (ms)');
    ylabel('Rate');
end

fprintf('\n--> [calc_psth.m]: spikechan  : %d\n', spikechan);
fprintf('--> [calc_psth.m]: syncchan   : %d\n', syncchan);
fprintf('--> [calc_psth.m]: duration_ms: %d\n', duration_ms);
fprintf('--> [calc_psth.m]: binwidth   : %d\n\n', binwidth);










