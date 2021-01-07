function h = revcor(h, t, ch)
% REVCOR - Compute reverse correlation between spike train and stimulus waveform
% REVCOR(H, T, CH) returns the revcor based on the histogram parameters in H
% and the time event data T/CH.
%

% get size statistics for stimulus, histogram and spikes
xsize = readUserFile(h.stimfile, 'size');
xlen = xsize(1);

dim = which_dim(h);
nbin = size(h.data,dim);

nspike = sum(ch == h.spikechan);

% Either of two methods is used: Cross correlation with PSTH for relatively
% long spike trains and short stimuli, and computing PESE then summing over spikes
% for short spike trains and long stimuli.  
% The factor of 1 is empirical.
%
% disp(nspike*nbin/(xlen*log2(xlen)))

if xlen*log2(xlen) < 1*nspike*nbin,
    
    disp('Using PSTH/xcorr method');
    tic
    
    % load stimulus and scale by reference voltage
    [x, Fs, ref] = readUserFile(h.stimfile, 'Fs', 'ref');
    if ~isnan(ref.Val),
        x = x/(10.^(.05*ref.Val));
    end    
    
    % histogram bin width must match stimulus sampling interval
    h.binwidth(dim) = 1000/Fs;    % assumes Fs in Hz
    lag0 = round(h.offset(dim)/h.binwidth(dim));
    lagend = lag0 + nbin;
    maxlag = max(abs([lag0 lagend]));
    
    % create and compute PSTH
    psth = set(h, 'Type', 'PST', 'Size', [xlen 1], 'BinWidth', 1000/Fs, 'Offset', 0);
    psth = hcompute(psth, t, ch);
    
    % repeat end of stimulus to handle case of continuous stimulus
    tsync = t(find(ch == h.syncchan));
    if abs((tsync(end)-tsync(1))/((psth.synccount-1)*h.binwidth(dim)) - xlen) < 3;
        x = [x(end-lagend+1:end); x];   % continuous stimulus
        psth.data = [zeros(lagend,1); psth.data];
    end
        
    % crosscorrelate stimulus with PSTH
    [xc, lags] = xcorr(psth.data, x, maxlag);
    h.data = xc(find(lags >= lag0 & lags < lagend));
    if dim == 2, h.data = h.data'; end

    h.synccount = psth.synccount;
    h.spikecount = psth.spikecount;
%    toc
    
else
    
    disp('Using PESE method');
%    tic

    pese = set(h, 'Type', 'PESE', 'Size', [nspike nbin], ...
               'Offsets', [1 h.offset(dim)]);
    
    pese = hcompute(pese, t, ch);
    if dim == 2,
        h.data = sum(pese.data);
    else
        h.data = sum(pese.data).';
    end
    h.binwidth(dim) = pese.binwidth(2);
    h.synccount = pese.synccount;
    h.spikecount = pese.spikecount;
%    toc
    
end

    
