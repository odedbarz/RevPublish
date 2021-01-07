function h = pese(h, t, ch)
% PESE - Get Pre-event stimulus ensemble (after Johannesma)
% PESE(H, T, CH) returns the PESE as a 2-D histogram
% Each row corresponds to one spike, each colum to a pre-spike time
% The Revcor is the mean of the PESE across spikes
%

% load stimulus and scale by reference voltage
[x, Fs, ref] = readUserFile(h.stimfile, 'Fs', 'ref');
x = x(:);
if ~isnan(ref.Val),
    x = x/(10.^(.05*ref.Val));
end
xlen = length(x);

% histogram bin width must match stimulus sampling interval
h.binwidth(2) = 1000/Fs;    % assumes Fs in Hz
nbin = size(h.data, 2);
lag0 = round(h.offset(2)/h.binwidth(2));
lagend = lag0 + nbin;

% get PS times and scale to multiples of bin width 
[tps, nsync] = pstimes(t, ch, h.syncchan);
kspike = find(ch == h.spikechan);
tps = tps(kspike);
tps = tps(find(~isnan(tps) & tps < xlen*h.binwidth(2)));
kps = floor(tps/h.binwidth(2)) + 1 + lagend;
    
% add stuff at beginning of stimulus to ease indexing
tsync = t(find(ch == h.syncchan));
if abs((tsync(end)-tsync(1))/((nsync-1)*h.binwidth(2)) - xlen) < 3;
    x = [x(end-lagend+1:end); x];   % continuous stimulus
else
    x = [zeros(lagend,1); x];       % bursts of sound
end
    
% use fancy indexing to get times preceding each spike
bins = [lag0:lagend-1];
kps = kps(:, ones(1,nbin)) - bins(ones(length(kps),1), :);
if size(kps,1) < size(h.data,1),
    h.data(1:size(kps,1),:) = x(kps);
else
    h.data = x(kps(1:size(h.data,1),:));
end
    
% additional information
h.binwidth(1) = 1;
h.offset(1) = 1;
h.synccount = nsync;
h.spikecount = length(kspike);
