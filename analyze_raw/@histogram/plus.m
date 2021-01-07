function hs = plus(h1, h2)
% HISTOGRAM/PLUS - Add two histograms 
% Sizes, types and bin widths must match
% One of the arguments can be a scalar
%
if isa(h1,'histogram') & isa(h2, 'histogram'),
    if ~strcmp(h1.type, h2.type)
        error('Histogram types do not match');
    end
    if any(h1.binwidth ~= h2.binwidth)
        error('Histogram binwidths do not match');
    end
    if any(size(h1.data) ~= size(h2.data))
        error('Histogram sizes do not match');
    end

    hs = h1;
    hs.data = h1.data + h2.data;
    hs.synccount = h1.synccount + h2.synccount;
    hs.spikecount = h1.spikecount + h2.spikecount;
    hs.periodcount = h1.periodcount + h2.periodcount;
    
elseif isa(h1, 'histogram'),
    hs = h1;
    hs.data = h1.data + h2;
elseif isa(h2, 'histogram'),
    hs = h2;
    hs.data = h1 + h2.data;
else
    error('At least one argument must be a histogram');
end
    
