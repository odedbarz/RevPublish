function hd = minus(h1, h2)
% HISTOGRAM/MINUS - Subtract two histograms 
% Sizes, types and bin widths must match
% Second argument can be a scalar
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

    hd = h1;
    hd.data = h1.data - h2.data;
    hd.synccount = h1.synccount - h2.synccount;
    hd.spikecount = h1.spikecount - h2.spikecount;
    hd.periodcount = h1.periodcount - h2.periodcount;
    
elseif isa(h1, 'histogram'),
    hd = h1;
    hd.data = h1.data - h2;
else
    error('First argument must be a histogram');
end
    