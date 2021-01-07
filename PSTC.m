function pst_st = PSTC(t, ch, duration_ms, binwidth, N, spikechan, syncchan, fignum)
%
%   pst_st = pstc(t, ch, duration_ms, binwidth, N, spikechan, syncchan, fignum)
%
%
% Description:
% Transforms spike trials to PSTH(s).
%

    % Make sure that t & ch are cells
    if ~iscell(t)   
        t = {t};
        ch = {ch};
        output_flag = 'array';
    else
        output_flag = 'cells';
    end


    if 4 > nargin || isempty(binwidth)
        binwidth  = 1;  % (ms)
    end

    if 5 > nargin || isempty(N)
        N  = 3;  % (ms)
    end

    if 7 > nargin || isempty(syncchan)
        syncchan  = 0;  % (ms)
    end

    if 9 > nargin || isempty(fignum)
        fignum = [];
    end
    
    % Set SPIKECHAN
    if 6 > nargin || isempty(spikechan)
        ch_list = cellfun(@(X) X(:)', ch, 'UniformOutput', false);
        ch_list = unique( [ch_list{:}] );
        ch_list = ch_list(ch_list > 0);     % avoid unneseccary syncchan (which are set to -1)
        ch_list = ch_list(ch_list ~= syncchan);
        assert(~isempty(ch_list), '--> ERROR at [PSTH.m]: there are no spikes in these measurements!!');
        assert(length(ch_list)==1, '--> ERROR at [PSTH.m]: there are TOO MANY spikes in these measurements -- you need to CHOOSE ONE!!');
        spikechan = ch_list(1);
    end
    
    % Run the smooth-PSTH for each cell
    [psth, bins, nsync_c, spikecount_c] = cellfun(@(T,CH)...
        pst_smooth_arrays(T, CH, duration_ms, binwidth, N, spikechan, syncchan, fignum), ...
        t, ch, 'UniformOutput', false);

    bins = bins{1};
    
    non_empty = cellfun(@(X) ~isempty(X), nsync_c);
    
    nsync = zeros(size(nsync_c));
    nsync(non_empty) = cellfun(@(X) X, nsync_c(non_empty));
    
    spikecount = zeros(size(spikecount_c));
    spikecount(non_empty) = cellfun(@(X) X, spikecount_c(non_empty));
    
    % If the input was a simple array, then return an array
    if strcmpi('array', output_flag)
        psth = psth{1}; 
    end

    % Set the output
    pst_st.type         = 'continuous psth';
    pst_st.duration_ms  = duration_ms;
    pst_st.binwidth     = binwidth;
    pst_st.fs           = 1/units.ms2sec(pst_st.binwidth);          % (Hz)    
    pst_st.N            = N;
    pst_st.spikechan    = spikechan;
    pst_st.syncchan     = syncchan;    
    pst_st.bins         = bins;
    pst_st.nsync        = nsync;
    pst_st.spikecount   = spikecount;
    pst_st.H            = [psth{:}];    % convert into a matrix
    
end



function [psth, bins, nsync, spikecount] = pst_smooth_arrays(t, ch, duration_ms, binwidth, N,...
    spikechan, syncchan, fignum)
%
% function [psth, bins, nsync, spikecount] = pst_smooth_arrays(t, ch, duration_ms, binwidth, N,...
%     spikechan, syncchan, fignum)
%
% This function calculates the smooth PSTH for T & CH array inputs.
%

    if isempty(t)
        psth = [];
        bins = [];
        nsync = [];
        spikecount = [];
        return;
    end

    histSize   = ceil(duration_ms/binwidth);
    bins       = (0 : (histSize-1)) * binwidth;
    tps        = spiketools.pstimes(t,ch);
    kspike     = find(ch == spikechan);
    nsync      = nnz(0==ch);            
    spikecount = nnz(kspike);        

    bins_plus_one = (0 : histSize) * binwidth;
    psth      = psth_smooth_one_trial(tps(kspike), bins_plus_one, N);
    psth      = psth(1:end-1);    % discard last, overflow bin

    % Normalize by the number of repeated trials
    if 0 ~= nsync
        psth = psth/nsync;
    end

    if isempty(fignum), return; end
    %%
    figure(99);
    clf;
    p1 = double(psth);
    plot(bins, p1)

    fprintf('\n--> [calc_psth.m]: spikechan  : %d\n', spikechan);
    fprintf('--> [calc_psth.m]: syncchan   : %d\n', syncchan);
    fprintf('--> [calc_psth.m]: duration_ms: %d\n', duration_ms);
    fprintf('--> [calc_psth.m]: binwidth   : %d\n\n', binwidth);

end


function psth = psth_smooth_one_trial(tps, bins, N)
%
%   function psth = psth_smooth_one_trial(tps, bins, N)
%
% Description:
% This function searches for N spikes (in TPS) at each time (BINS). Finding
% N spikes at a given time defines the symmetric window size for this time.
%

    psth = [];
    if isempty(tps), return; end
    if 3 > nargin, N = 3; end

    bins = bins(:);
    len_psth = length(bins);
    tps = sort(tps, 'ascend');
    win_size_ms = zeros(len_psth, 1);
    %psth = nan(len_psth,1);
    
    for nn = 1:len_psth
        t0 = bins(nn);  % the n'th time bin
        
        % Option #1: Sort (TOO slow!)
        %dist_abs = sort(abs(tps - t0), 'ascend');
        %win_size_ms(nn) = dist_abs(N);        % (ms)
        
        % Option #2: avoid the sort
        dist_abs = abs(tps - t0);
        for mm = 1:N
            [w_ms, idx_k] = min(dist_abs);
            if idx_k > len_psth, break; end            
            dist_abs(idx_k) = Inf;
        end
        win_size_ms(nn) = w_ms;        % (ms)
        
    end
    
    % Divide by the window size at each time location 
    psth = N./units.ms2sec(win_size_ms);   % (1/ms) --> (Hz=1/sec)
    

end







