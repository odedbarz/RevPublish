function hq = rdivide(h, s)
% HISTOGRAM/RDIVIDE - Divide a histogram by a scalar 
% Second argument must be a scalar
%
if isa(h,'histogram') & isa(s, 'histogram'),
    error('Can''t divide two histograms')    
elseif isa(h, 'histogram'),
    hq = h;
    hq.data = h.data ./ s;
else
    error('First argument must be a histogram');
end
    