function psth_st = PSTW(S, duration_ms, binwidth, win_size_ms, spikechan, syncchan, fignum)
%
%   pst_st = pstw(t, ch, duration_ms, binwidth, spikechan, syncchan, fignum)
%
% Description:
% Transforms spike trials to PSTH(s).
%
% h holds the spike counts at each bin

% win_type = 'hanning';   % (str) window type for the smppthing
% win_size_ms = 25;       % (ms) window size

t = S.t;
ch= S.ch;

if 3 > nargin || isempty(binwidth)
    binwidth  = 1;  % (ms)
end

if 4 > nargin || isempty(win_size_ms)
    win_size_ms = 25;       % (ms) window size
end

if 6 > nargin || isempty(syncchan)
    syncchan  = 0;  % (ms)
end

if 5 > nargin || isempty(spikechan)
    %spikechan  = 1;  % (ms)    
    ch_list = unique( [ch{1}] );
    ch_list = ch_list(ch_list ~= syncchan);
    assert(~isempty(ch_list), '--> ERROR at [PSTH.m]: there are no spikes in these measurements!!');
    assert(length(ch_list)==1, '--> ERROR at [PSTH.m]: there are TOO MANY spikes in these measurements -- you need to CHOOSE ONE!!');
end

% if 8 > nargin || isempty(use_hcompute)
%     use_hcompute = true;
% end
% When performing windowing, do it over each PSTH and THEN do the
% averaging. That's why use_hcompute is always FALSE
use_hcompute = false;

if 7 > nargin || isempty(fignum)
    fignum = [];
end



%%
% When performing windowing, do it over each PSTH and THEN do the averaging
psth_st = PSTH(S, duration_ms, binwidth, spikechan, syncchan, use_hcompute, fignum);

% Define the window
psth_st.type       = 'windowed psth';
psth_st.win_size_ms= win_size_ms;       % (ms) window size
psth_st.n_win      = fix(psth_st.fs * 1e-3*psth_st.win_size_ms);     % (smp) window size; # of samples
% psth_st.win_fun    = @(n) hanning(n)/sum(hanning(n));
psth_st.win_fun    = @(n) gausswin(n)/sum(gausswin(n));

% Smooth PSTH using HANNING window
% psth_st.win      = hanning(psth_st.n_win)/sum(hanning(psth_st.n_win));
%psth_st.win      = psth_st.win_fun(psth_st.n_win)/sum(psth_st.win_fun(psth_st.n_win));
psth_st.win      = psth_st.win_fun(psth_st.n_win);
psth_st.Hbin     = psth_st.H;                         % response (spikes/sec)

%% Perform the windowing
% FILTFILT to avoids the lag
%psth_st.H = filtfilt(psth_st.win, 1, psth_st.Hbin);    % response (spikes/sec)

if isfield(psth_st, 'H_trials_win')
    psth_st.H_trials_win = cell(size(psth_st.H_trials));
    for q = 1:length(psth_st.H_trials_win(:))
        if 1 < psth_st.n_win
            psth_st.H_trials_win{q} = filtfilt(psth_st.win, 1, psth_st.H_trials{q});    % response (spikes/sec)
        else
            psth_st.H_trials_win{q} = psth_st.H_trials{q};
        end
    end
    
    % Get the STD around the average
    psth_st.H_SE = cellfun(@(AVG,X) X-AVG, psth_st.H, psth_st.H_trials_win, 'UniformOutput', false);
end


%psth_st.H = filtfilt(psth_st.win, 1, psth_st.H_trials);    % response (spikes/sec)
% Average over the windowed trials PSTHs of each case
% psth_st.H = cellfun(@(X) mean(X,2), psth_st.Hbin, 'UniformOutput', false);

% conv(\Sigma) <--> \Sigma(conv)
psth_st.H = arrayfun(@(I) filtfilt(psth_st.win, 1, psth_st.Hbin(:,I)), 1:size(psth_st.Hbin,2),...
    'UniformOutput', false);
psth_st.H = [psth_st.H{:}];     % Cell --> matrix


if isempty(fignum), return; end
%%
figure(fignum + 10);
clf;
p1 = double(psth_st{1});
plot(bins, p1)

fprintf('\n--> [calc_psth.m]: spikechan  : %d\n', spikechan);
fprintf('--> [calc_psth.m]: syncchan   : %d\n', syncchan);
fprintf('--> [calc_psth.m]: duration_ms: %d\n', duration_ms);
fprintf('--> [calc_psth.m]: binwidth   : %d\n\n', binwidth);










