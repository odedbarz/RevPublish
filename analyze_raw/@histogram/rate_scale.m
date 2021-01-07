function sf = rate_scale(h)
% RATE_SCALE - Get scale factor to express histogram in units of discharge
% rate (spikes/sec).  Most meaningful for 1-D histograms.
%
if num_dim(h.type) == 1,   % 1-D histogram
    dim = which_dim(h);
    bwdim = dim - ndims(h.data) + 2;
end

switch h.type
    case {'PST', 'Latency'}
        sf = 1000./(h.binwidth(bwdim)*h.synccount);
    case {'Interval', 'Autocorrelation', 'Crosscorrelation'}
        sf = 1000./(h.binwidth(bwdim)*h.spikecount(1));
    case {'Period'}
        sf = 1000*size(h.data,dim)./(h.period*h.periodcount);
    case 'Revcor'   % divide by total stimulation time
        stim_size = readUserFile(h.stimfile, 'size');
        sf = 1000./(h.synccount * stim_size(1) * h.binwidth(bwdim));
    case {'PST Raster'}
        sf = 1000./h.binwidth(2);
    case {'Interval Raster', 'Autocorrelation Raster', 'Crosscorrelation Raster'}
        sf = 1000*h.synccount./(h.binwidth(2)*h.spikecount(1));
    case {'Period Raster'}
        sf = 1000*h.synccount*size(h.data,ndims(h.data))./(h.period*h.periodcount);
    case {'PST-Interval', 'PST-Autocorrelation', 'PST-Crosscorrelation', 'PST-Period'}
        sf = 1e6./(h.synccount*prod(h.binwidth));
    case {'Period-Interval'}
        sf = 1e6./(h.periodcount*prod(h.binwidth));
    case {'Joint Interval' }
        sf = 1e6./(h.spikecount*prod(h.binwidth));
    otherwise
        sf = 1;
end

        
    
