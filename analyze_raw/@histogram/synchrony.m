function s = synchrony(h, scaleopt)
% SYNCHRONY - Compute the synchronized rate from period histogram
% S=SYNCHRONY(H) returns the comlex synchronized rate (actually synchronized spike count)
% fron the histogram H.  The magnitude of the comlex synchronized rate is the
% conventional synchronized rate, while its phase is the response phase.
%
% S= SYNCHRONY(H, 'index') returns the complex synchronization index , which is
% the synchronized count divided by the total spike count
%
% If H is a 1-D period histogram, S is a scalar.
% If H is a 2-D histogram one dimension of which is period (e.g a Period Raster or
% PST-Period histogram or a Period neurogram), S is a 1-D histogram containing 
% the synchronized rates for each value of the other dimension
%
if isempty(strfind(h.type, 'Period')),
    error('At least one dimension of histogram must be period')
end

switch h.type
case 'Period'
    dim = which_dim(h);
case 'Period Raster'
    dim = ndims(h.data);
    h.type = 'Synchrony Raster';
case 'PST-Period'
    dim = ndims(h.data)-1;
    h.type = 'PST-Synchrony';
case 'Period-Interval'
    dim = ndims(h.data);
    h.type = 'Synchrony-Interval';
otherwise
    error('Unsupported histogram type');
end

t = get(h, 'BinTimes', dim);
%   changed to make it compatible with normalized scale of PH
%   exp_vec = exp(-j*2*pi*t/h.period);  
exp_vec = exp(-1j*2*pi*t);
s = h;

if dim == 1
    s.data = exp_vec * h.data;
else
    s.data = full(h.data * exp_vec.');
end

if nargin >= 2 && strcmp(scaleopt, 'index')
    s.data = s.data ./ sum(h.data, dim);
end

if max(size(s.data)) == 1  % scalar
    s = s.data;
else
    s.binwidth(dim-ndims(h)+2) = h.period;
end
