function [delays, tps] = shuffled_intervals(t, ch, syncChan, spikeChan)
% [delays t] = shufled_intervals(t, ch, syncChan, spikeChan)
% Counts all delays between spike pairs in non-identical stimulus
% presentations, for the SAC/XAC analysis a la Joris, 2004.  Each
% trial, defined by the events on syncChan, is compared with all other
% tials excluding itself.
%
% There are two ways to compute XAC:
% 1. If the two spike trains to be crossed are from different stimuli, use
% cell vectors: t = {t1, t2}, ch = {ch1, ch2}.
% syncChan and spikeChan can be 1- or 2-element vectors, specifying the
% channels used in ch1 and ch2.
% 2. If the two spike trains are from the same stimulus, simultaneous
% recording on different channels, t and ch should be vectors containing
% all channel information.  In this case, syncChan must be a single element
% and spikeChan must be a 2-element vector denoting the two channels of
% spikes to be crossed.
%
% delays is returned unsorted
% tps (optional) is a vector of post-stimulus times for each interval

import spiketools.*

if isempty(ch)
    delays = [];
    if nargout > 1, tps = []; end
    return
end

if nargin < 3, syncChan=0; end

if ~iscell(t) && ~iscell(ch)
    
    if nargin < 4, spikeChan = 1 ; end
    
    if numel(spikeChan) == 1            %Single channel, SAC
        
        t = pstimes(t,ch,syncChan);
        %sort data into relevant vectors
        tSpike = t(ch==spikeChan);
        trialNum = cumsum(ch == syncChan);
        trialNum = trialNum(ch == spikeChan);

        %preallocate for speed, figure out how many intervals we have
        np = 0;
        for trial=1:trialNum(end)
            np = np+sum(trialNum==trial)*sum(trialNum>trial);
        end
        delays = zeros(np,1);
        if nargout > 1, tps = zeros(size(delays)); end
        delayIndex = 1;

        %calculate delays to all spikes in trials with higher trial number
        %not restricting to higher number would simply give symmetry as X->Y
        %has the the same but opposite delays as Y->X
        for spike=1:length(tSpike)
            kgt = trialNum > trialNum(spike);
            %figure out how many spikepairs we're dealing with
            numPairs = sum(kgt);
            if numPairs>0
                delays(delayIndex:delayIndex+numPairs-1) = tSpike(kgt) - tSpike(spike);
                if nargout > 1, tps(delayIndex:delayIndex+numPairs-1) = tSpike(kgt); end
                delayIndex = delayIndex + numPairs;
            end        
        end
        
    else     %multiple-channel, XAC

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
    
        
    end
    
elseif iscell(t) && iscell(ch)      %Two-stimulus XAC
    
    if nargin < 4, spikeChan=[1 1]; end
    
    t1 = pstimes(t{1},ch{1},syncChan);
    t2 = pstimes(t{2},ch{2},syncChan);
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
    
else
    error('t and ch must either be both vectors or both cell arrays!');
end
