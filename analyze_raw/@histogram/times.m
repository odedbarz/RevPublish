function hp = times(h1, h2)
% HISTOGRAM/TIMES - Multiply a histogram by a scalar 
% One of the arguments must be a scalar
%
if isa(h1,'histogram') & isa(h2, 'histogram'),
    error('Can''t multiply two histograms')    
elseif isa(h1, 'histogram'),
    hp = h1;
    hp.data = h1.data .* h2;
elseif isa(h2, 'histogram'),
    hp = h2;
    hp.data = h1 .* h2.data;
else
    error('At least one argument must be a histogram');
end
    