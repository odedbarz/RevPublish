function [ti, kspike] = intervals(t, ch, spikechan, order)
% INTERVALS - Compute interspike intervals from event times
% INTERVALS(T, CH, SPIKECHAN) returns the first-order intervals 
% from the event times T and channels CH for spikes in channel SPIKECHAN (default 1)
%
% INTERVALS(T, CH, SPIKECHAN, ORDER) returns the intervals of the specified ORDER
%
% [TI, KSPIKE] = INTERVALS(...) returns the intervals TI as well as a NINERTVAL X 2 matrix 
% KSPIKE of indices to the first and second spike forming each interval, 
% i.e. the intervals are T(KSPIKE(:,2)) - T(KSPIKE(:,1))
%
%

if nargin < 3, spikechan = 1; end

if nargin < 4
    order = 1; 
elseif order <= 0
    error('Interval order must be > 0');
end

kspike = find(ch == spikechan);
ts = t(kspike);

ti = ts(order+1:end) - ts(1:end-order);

%{
if nargout > 1
    kspike = [kspike(1:end-order) kspike(order+1:end)];
end
%}

