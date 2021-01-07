function [delays, tps] = shuffled_cross_intervals(t, ch, syncChan, spikeChan)
% delays = shuffled_cross_intervals(t, ch, syncChan, spikeChan)
% like shuffled_intervals, except ignore pairings between identical trials
% spikeChan must contain two elements, referring to the channel numbers for
% the spikes to be cross-correlated
% This function is used for simultaneous two-unit recording where the
% response of either unit to the same trial is coupled by network
% connections.  Those effects are discarded by ignoring trial pairs between
% identical trials, which is normally not done for the SXC.
%
% delays is returned unsorted
% tps (optional) is a vecotr of post-stimulus times for each interval

%force column dimnesion
t = t(:);
ch = ch(:);

if ~iscell(t)
    %input has two channels of spikes to cross
    if nargin < 3, syncChan = 0; end
    if nargin < 4, spikeChan = [1 2]; end

    if numel(spikeChan) ~= 2, 
        error('Need either two spike channels or a cell array of two spike variables for shuffled cross intervals');
    end
    
    t = pstimes(t,ch,syncChan);

    %sort data into relevant vectors
    tSpike1 = t(ch==spikeChan(1));
    tSpike2 = t(ch==spikeChan(2));
    trialNum = cumsum(ch == syncChan);
    trialNum1 = trialNum(ch==spikeChan(1));
    trialNum2 = trialNum(ch==spikeChan(2));

    %preallocate for speed, figure out how many intervals we have
    np = 0;
    for trial=1:sum(ch==syncChan)
        np = np+sum(trialNum1==trial)*sum(trialNum2~=trial);
    end
    delays = zeros(np,1);
    if nargout > 1, tps = zeros(size(delays)); end
    delayIndex = 1;

    %calculate delays to all spikes in non-identical trials
    for spike1=1:length(tSpike1)
        kne = trialNum2 ~= trialNum1(spike1);
        %figure out how many spikepairs we're dealing with
        numPairs = sum(kne);
        if numPairs>0
            delays(delayIndex:delayIndex+numPairs-1) = tSpike2(kne) - tSpike1(spike1);
            if nargout > 1, tps(delayIndex:delayIndex+numPairs-1) = tSpike2(kne); end
            delayIndex = delayIndex + numPairs;
        end        
    end
    
else
    if nargin < 3, syncChan=0; end
    if nargin < 4, spikeChan=1; end

    t1 = pstimes(t{1},ch{1},syncChan(1));
    t2 = pstimes(t{2},ch{2},syncChan(end));
    %sort data into relevant vectors
    tSpike1 = t1(ch{1}==spikeChan(1));
    tSpike2 = t2(ch{2}==spikeChan(end));

    %preallocate for speed
    delays = zeros(size(tSpike1).*size(tSpike2));
    if nargout > 1, tps = zeros(size(delays)); end

    %for each spike in train1, find delay to all spikes in train2.
    for spike = 1:length(tSpike1)
        delays((spike-1)*length(tSpike2)+1:spike*length(tSpike2)) = tSpike2 - tSpike1(spike);
        if nargout > 1, % Return PS times for train2
            tps((spike-1)*length(tSpike2)+1:spike*length(tSpike2)) = tSpike2;
        end
    end
end