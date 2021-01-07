function [ti, kspike] = all_order_interv(t, ch, spike_chan, maxinterv)
% ALL_ORDER_INTERV - Compute all-order interspike intervals from event times
% ALL_ORDER_INTERV(T, CH, SPIKECHAN, MAXINTERV) computes the all-order intervals 
% from the event times T and channels CH for spikes in channel SPIKECHAN (default 1)
% up to the largest interval MAXINTERV (default Inf)
%
%  [TI, KSPIKE] = ALL_ORDER_INTERV(T, CH, SPIKECHAN, MAXINTERV) returns
%  indices to the spikes in each interval as well as the intervals TI
%  KSPIKE is an NINTERVALS X 2 matrix containing the index of the first and
%  second spike forming each interval, i.e. the intervals are
%  T(KSPIKE(:,2))-T(KSPIKE(:,1))
%
if nargin < 3, spike_chan = 1; end

% select spikes only
ks = find(ch == spike_chan);
if isempty(ks),
    ti = [];
    kspike = [];
    return;
end
ts = t(ks);
ts = ts(:);
ns = length(ts);

% For relatively short spike trains and long maximum intervals,
% we use an algorithm with no loop which creates an [nspike X nspike] matrix.
% Otherwise, we use an algorithm wich loops over the interval order
% The 30 factor is based on empirical speed tests
%
if nargin < 4, maxinterv = ts(end) - ts(1); end
mean_interv = (ts(end) - ts(1))/ns;
mean_order = maxinterv/mean_interv;
% disp(ns)
% disp(mean_order)
% disp(mean_interv)

if ns < 3 * mean_order,     % use full matrix
%    disp('using full matrix method')
    ti = ts(:,ones(1,ns));
    ti = ti - ti';
    klt = find(ti > 0 & ti < maxinterv);
    ti = ti(klt);
    
    if nargout > 1,
        [irow,icol] = ind2sub([ns ns], klt);
        kspike = [ks(icol) ks(irow)];
    end
    
else,           % loop over interval order
%    disp('using loop method')
    ti = [];
    kspike = [];

    for order = 1:ns-1;
        to = ts(order+1:end) - ts(1:end-order);
        if all(to >= maxinterv), break; end
        ti = [ti; to];
        if nargout > 1,
            kspike = [kspike; [ks(1:end-order) ks(order+1:end)]];
        end
    end
end



