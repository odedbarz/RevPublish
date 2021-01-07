function [tps, ch_new, t_new] = get_spike_times(t, ch, spikechan, syncchan)
%
%   function tps = get_spike_times(t, ch, spikechan, [syncchan])
%
%
% Description: 
% Extract the spike times.
%

if ~exist('syncchan', 'var')
    syncchan = 0;
end

% In this case, ALL inputs are double vectors and the output will be too
if ~iscell(t)
    double_flag = true;
    t = {t};
    ch = {ch};
else
    double_flag = false;
end


%%
% Extract the selected spikechan
spikes = cellfun(@(CH) CH==spikechan, ch, 'UniformOutput', false);
chans  = cellfun(@(CH) CH==syncchan, ch, 'UniformOutput', false);

% Keep the desired spike-channel and sync-channel
t_new  = cellfun(@(T,IDX1,IDX2) T(IDX1 | IDX2), t, spikes, chans, 'UniformOutput', false);
ch_new = cellfun(@(CH,IDX1,IDX2) CH(IDX1 | IDX2), ch, spikes, chans, 'UniformOutput', false);

% Spike times
[tps, ~] = cellfun(@(T,CH) spiketools.pstimes(T,CH), t_new, ch_new, 'UniformOutput', false);


%% Output
if true == double_flag
    t_new = t_new{1};
    ch_new = ch_new{1};
    tps = tps{1};
end











