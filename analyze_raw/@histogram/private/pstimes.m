function [tps, nstim] = pstimes(t, ch, syncchan)
% PSTIMES - Compute poststimulus times of spike data
% Usage: [tps, nsync] = pstimes(event_times, event_channels, syncchan)
%       event_times     vector of event times in milliseconds
%       event_channels  vector of event channels
%       syncchan    synch-pulse channel (default 0)
%	    tps		    post-stimulus times of ALL the events in msec
%
if nargin < 3, syncchan = 0; end
t = t(:);

% find sync times
ind_sync = find(ch==syncchan); 
nstim = length(ind_sync);

% treat special case separately
if nstim == 0,
    tps = zeros(size(t));
    tps(:) = NaN;
    return;
end

% create an array of delta functions at the sync times, 
% each weighted by increment in PST
ts = zeros(size(t));
ts(ind_sync) = diff([0; t(ind_sync)]);

% subtract step function at sync times from original
tps = t - cumsum(ts);

% PST is undefined before first sync pulse
tps(1:ind_sync(1)-1) = NaN;
