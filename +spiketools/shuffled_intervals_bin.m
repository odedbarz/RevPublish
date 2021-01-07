function binnedData = shuffled_intervals_bin(t, ch, binWidth, numBins, binOffset, syncChan, spikeChan)
% delays = shufled_intervals(t, ch, syncChan, spikeChan)
% Counts all delays between spike pairs in non-identical stimulus
% presentations, for the SAC analysis a la Joris, 2004.  Each
% trial, defined by the events on syncChan, is compared with all other
% tials excluding itself.
%
% if syncChan and spikeChan are not specified, 0 and 1 are used.
%
% To compute XAC delays for two stimuli, use
% shuffled_intervals({t1 t2}, {ch1 ch2}, syncChan, spikeChan)
% syncChan and spikeChan can each be one or two element vectors.  If one
% element, the same channel is used for each t,ch.
%
% Results are binned according to binWidth, numBins, binOffset.  If 

if nargin<6
    syncChan = 0;
    spikeChan = 1;
end

if nargin<5
    binOffset = [0 0];
end

dWidth = binWidth(1);
numD = numBins(1);
dOffset = binOffset(1);

if numel(numBins)==2 && numel(binWidth)==2
    if numel(binOffset)==1
        binOffset(2) = 0;
    end
    tWidth = binWidth(2);
    numT = numBins(2);
    tOffset = binOffset(2);
else  %use a single time bin
    if iscell(t)
        tWidth = max(t{1}(end),t{2}(end));
    else
        tWidth = t(end);
    end
    numT = 1;
    tOffset = 0;
end

binnedData = zeros(numD,numT);
    
if ~iscell(t) %SAC

    t = pstimes(t,ch,syncChan);
    %sort data into relevant vectors
    tSpike = t(ch==spikeChan);
    trialNum = cumsum(ch == syncChan);
    trialNum = trialNum(ch==spikeChan);
    
    N = sum(ch==syncChan);
    trialNorm = N*(N-1)/2;

    if ~isempty(trialNum)
        for trial=1:trialNum(end)-1
            otherSpikes = tSpike(trialNum > trial);  %~= would just give symmetry
            theseSpikes = tSpike(trialNum == trial);
            for spike=1:sum(trialNum==trial)
                delays = otherSpikes-theseSpikes(spike);
%                times = max(otherSpikes, theseSpikes(spike));
                
                % find t & d indices
                dI = ceil((delays-dOffset)/dWidth);
%                tI = ceil((times-tOffset)/tWidth);
                tI = ceil(theseSpikes(spike)/tWidth);
                %as long as the max delay << binWidth, these are more or
                %less equivalent

                % restrict indices to matrix dimensions
                ok = (dI >= 1 & dI <= numD & ...
                      tI >= 1 & tI <= numT);

                % use the Jay Delosge trick
                data = sparse(dI(ok),tI,1, numD,numT);
                %if we use times to compute tI, we must use tI(ok)
                
                binnedData = binnedData + data; 
            end
        end
    end
else      % XAC
    t1 = psTimes(t{1},ch{1},syncChan(1));
    t2 = psTimes(t{2},ch{2},syncChan(end));

    %sort data into relevant vectors
    tSpike1 = t1(ch{1}==spikeChan(1));
    tSpike2 = t2(ch{2}==spikeChan(end));

    N1 = sum(ch{1}==syncChan(1));
    N2 = sum(ch{2}==syncChan(end));
    trialNorm = N1*N2;

    
    %for each spike in train1, find delay to all spikes in train2.
    for spike = 1:length(tSpike1)
        %binnedDelays = binnedDelays + histc(tSpike2-tSpike1(spike),binEdges,1);
        
        delays = tSpike2-tSpike1(spike);
%        times = max(tSpike2, tSpike1(spike));

        % find t & d indices
        dI = ceil((delays-dOffset)/dWidth);
%        tI = ceil((times-tOffset)/tWidth);
        tI = ceil(tSpike1(spike)/tWidth);

        % restrict indices to matrix dimensions
        ok = (dI >= 1 & dI <= numD & ...
              tI >= 1 & tI <= numT);

        % use the Jay Delosge trick
        data = sparse(dI(ok),tI,1, numD,numT);

        %binnedDelays = binnedDelays + histc(delays,binEdges,1);
        binnedData = binnedData + data;
    end
end

%normalize
for tBin = 1:numT
    binnedData(:,tBin) = binnedData(:,tBin) / (trialNorm * sum(binnedData(:,tBin)) / tWidth);
end