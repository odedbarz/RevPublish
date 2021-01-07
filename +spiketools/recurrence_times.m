function [ti, kspike] = recurrence_times(t, ch, spike_chan, direction)
% RECURRENCE_TIMES - Compute recurrence times between two simultaneous spikes trains
% RECURRENCE_TIMES(T, CH, SPIKECHAN, DIRECTION) returns the recurrence times between
% the spikes in SPIKECHAN(1) and those in SPIKECHAN(2).  DIRECTION is either 'forward' 
% 'backward' or 'both' (the default).  Forward recurrence times
% are returned with a positive sign, backward recurrence times with a negative sign.
% If SPIKECHAN is a scalar, single channel first-order intervals are
% returned.
%
% [TI, KSPIKE] = RECURRENCE_TIMES(...) returns the intervals TI as well as
% a NINERTVAL X 2 matrix of indices to the first and second spike forming
% each interval, i.e. the recurrence times are T(KSPIKE(:,2)) - T(KSPIKE(:,1))
%
if nargin < 4, direction = 'both'; end

if nargin < 3, spike_chan = [1 2]; 
elseif length(spike_chan) < 2, spike_chan = [spike_chan spike_chan]; 
end

t = t(:);   % make column vectors
ch = ch(:);

% select spike times
ks = find(ch == spike_chan(2));
dmax = max(diff(ks))-1;
kf = find(ch == spike_chan(1)); % forward indices
kb = kf;    % backward indices
kspike = [];

switch lower(direction)
    case {'forward', 'both'}
        for d = 1:dmax,    % loop over delta
            kf(find(kf+d > ks(end))) = [];   % trim end to avoid indexing errors
            keq = find(ch(kf+d) == spike_chan(2));
            kgood = kf(keq);
            kspike = [kspike; [kgood kgood+d]]; 
            kf(keq) = []; % remove done items
        end
end

switch lower(direction)
    case {'backward', 'both'}
         for d = 1:dmax,
            kb(find(kb-d < ks(1))) = [];   % trim beginning to avoid indexing errors
            keq = find(ch(kb-d) == spike_chan(2));
            kgood = kb(keq);
            kspike = [kspike; [kgood kgood-d]];
            kb(keq) = [];
        end
end

ti = diff(t(kspike),1,2);