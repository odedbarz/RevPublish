function [tm, chm] = merge_trials(t, ch, sync_chan)
% MERGE_TRIALS _ Merge all stimulus trials into one where spikes are
% ordered by peri-stimulus times
% Usage: [tm, chm] = merge_trials(t, ch, sync_chan)
% t, ch:  Original spike train
% sync_chan: Optional sync channel (default 0)
% tm, chm: merged spike train
%
if nargin < 3, sync_chan = 0; end

[tm, nsync] = pstimes(t, ch, sync_chan);
chm = ch;

% remove spikes before first sync pulse
knan = find(isnan(tm));
tm(knan) = [];
chm(knan) = [];

if nsync <= 1; return; end

% sort spikes
[tm, ks] = sort(tm);
chm = chm(ks);

% remove extra sync pulses
tm(1:nsync-1)=[];
chm(1:nsync-1)=[];




