function [ti, kspike] = intervals(t, ch, spike_chan, order)
% INTERVALS - Compute interspike intervals from event times
% INTERVALS(T, CH, SPIKECHAN) returns the first-order intervals 
% from the event times T and channels CH for spikes in channel SPIKECHAN (default 1)
%
% INTERVALS(T, CH, SPIKECHAN, ORDER) returns the intervals of the specified ORDER
%
% [TI, KSPIKE] = INTERVALS(...) returns the intervals TI as well as a vector
% of index to the spike times in the original array T
%
if nargin < 4, 
    order = 1; 
elseif order <= 0
    error('Interval order must be > 0');
end
if nargin < 3, spike_chan = 1; end

kspike = find(ch == spike_chan);
ts = t(kspike);
tbad = NaN;

ti = [tbad(ones(order,1)); ts(order+1:end) - ts(1:end-order)];


