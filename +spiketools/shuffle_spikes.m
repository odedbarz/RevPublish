function [ts, chs] = shuffle_spikes(t, ch, syncchan, spikechan)
% SHUFFLE_SPIKES(T, CH, SYNCCHAN, SPIKECHAN) creates a virtual array of
% spike times for computation of shuffed autocorrelations a la Joris (2004).
% If the input array has N trials, the output array will have (N X N-1)/2 pseudo trials
% consisting of combined spikes from every distinct pair of trials in the
% input.  Output channels are SPIKECHAN and SPIKECHAN+1.
%

import spiketools.*

if nargin < 3, syncchan = 0; end
if nargin < 4, spikechan = 1; end

% keep only necessary data
keep = find(ch == syncchan | ch == spikechan);
t = t(keep);
ch = ch(keep);

% find trial onsets and PS times
ksync = find(ch == syncchan);
nsync = length(ksync);
tdiff = max(diff(t(ksync)));    % pseudo time increment between trials
trialno = cumsum(ch == syncchan);
tps = pstimes(t, ch, syncchan);

% add bogus trial for ease of indexing
ksync = [ksync; length(t)+1];

% first loop to identify total number of spikes
ns = 0;
for kdiff = 1:nsync-1,  % loop over trial offsets
    n1 = ksync(nsync+1-kdiff)-ksync(1);
    n2 = ksync(end)-ksync(kdiff+1) - nsync + kdiff;
    ns = ns + n1 + n2;
end

% preallocate arrays for speed
ts = zeros(ns,1);
chs = zeros(ns,1);
toffset = 0; % time offset
noff = 0;   % offset in array

for kdiff = 1:nsync-1,  % loop over trial offsets
    k1 = ksync(1):ksync(nsync+1-kdiff)-1;
    k2 = ksync(kdiff+1):ksync(end)-1;
    k2 = k2(find(ch(k2)==spikechan));   % keep only spikes
    nadd = length(k1)+length(k2);
    t1 = tps(k1) + (trialno(k1)-1)*tdiff + toffset;
    t2 = tps(k2) + (trialno(k2)-kdiff-1)*tdiff + toffset;
    [ts(noff+1:noff+nadd), ks] = sort([t1; t2]);
    chadd = [ch(k1); ch(k2)+1];
    chs(noff+1:noff+nadd) = chadd(ks);
    toffset = toffset + tdiff * (nsync-kdiff);   
    noff = noff + nadd;
end


