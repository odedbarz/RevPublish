function [ti, kspike] = cross_intervals(t, ch, spikechan, maxinterv)
% CROSS_INTERVALS - Find all order intervals between two spike trains
% CROSS_INTERVALS(T, CH, SPIKECHAN, MAXINTERV) returns the intervals between
% the spikes in SPIKECHAN(1) and SPIKECHAN(2) up to a maximum MAXINTERV (default Inf)
% If SPIKECHAN is a scalar, single channel all-order intervals are returned
%
% [TI, KSPIKE] = CROSS_INTERVALS(T, CH, SPIKECHAN, MAXINTERV) returns
% indices to the spikes in each interval as well as the intervals TI
% KSPIKE is an NINTERVALS X 2 matrix containing the index of the first and
% spikes spike forming each interval, i.e. the cross intervals are
% T(KSPIKE(:,2)) - T(KSPIKE(:,1))
%
if nargin < 3, spikechan = [1 2]; end
if nargin < 4, maxinterv = t(end) - t(1); end
if length(spikechan) < 2; spikechan = [spikechan spikechan]; end

t = t(:);   % make column vectors
ch = ch(:);

% For relatively short spike trains and long maximum intervals,
% we use an algorithm with no loop which creates an [nspike X nref] matrix.
% Otherwise, we use an algorithm wich loops over the interval order.
% The 30 factor is based on empirical speed tests
%
mean_interv = (t(end) - t(1))/length(t);
mean_order = maxinterv/mean_interv;

n1 = sum(ch == spikechan(1));
n2 = sum(ch == spikechan(2));

if n1 * n2 < 18 * mean_order * (n1 + n2),  % use full matrix method
%   tic

%    display('using full matrix method')
    % get spike times in each channel
    k1 = find(ch == spikechan(1));
    k2 = find(ch == spikechan(2));
    t1 = t(k1);   % column vector
    t2 = t(k2)';  % row vector
   
    % subtract the spike times pairwise to get the intervals
    ti = t2(ones(n1,1), :) - t1(:, ones(1,n2));

    % eliminate intervals out of range 
    klt = find(abs(ti) < maxinterv);
    ti = ti(klt);

    % compute indices
    if nargout > 1,
        [irow,icol] = ind2sub([n1 n2], klt);
        kspike = [k1(irow) k2(icol)];  % [ninterv X 2] matrix
    end
    
%    toc
else        % use loop method
%    tic
%   display('using loop method')
    
    % find spike times from either channel
    ks = find(ch == spikechan(2));
  
    if spikechan(1) == spikechan(2),    %include zero-order intervals for ACH
        kf = ks;
    else
        kf = find(ch == spikechan(1));  % forward indices
        kb = kf;
    end
    kspike = [];
    
%   do forward intervals first
    for d = 1:ks(end)-kf(1),     % loop over increment
        kf(find(kf+d > ks(end))) = [];   % trim end to avoid indexing errors
        kf(find(t(kf+d)-t(kf) > maxinterv)) = [];
        if isempty(kf), break; end
        keq = find(ch(kf+d) == spikechan(2));
        kgood = kf(keq);
        kspike = [kspike; [kgood kgood+d]]; 
    end
    
    % now, backward intervals
    if spikechan(1) == spikechan(2),
        kspike = [kspike; [ks ks]; fliplr(kspike)];
    else
        for d = 1:kb(end)-ks(1),     % loop over increments
            kb(find(kb-d < ks(1))) = [];   % trim beginning to avoid indexing errors
            kb(find(t(kb)-t(kb-d) > maxinterv)) = [];
            if isempty(kb), break; end
            keq = find(ch(kb-d) == spikechan(2));
            kgood = kb(keq);
            kspike = [kspike; [kgood kgood-d]];
        end
    end
    
    ti = diff(t(kspike),1,2);
%    toc
end

