function hs = sum(h, dim)
% HISTOGRAM/SUM - Sum a histogram
% For 1-D histograms, SUM(H) returns a scalar equal to the sum of the data in H
% For 2-D histogram, SUM(H) returns a 1-D histogram containing the sum of
% each column of the data in H
% SUM(H, DIM) returns the sum along dimension DIM
%
if nargin == 1, dim = which_dim(h); end

if min(size(h.data)) == 1,  % 1-D histogram
    hs = sum(h.data, dim);
else    % 2-D histogram
    hs = h;
    hs.data = full(sum(h.data, dim));
    % collapse binwidths
    hs.binwidth(dim) = h.binwidth(dim)*size(h.data,dim);
    
    % guess type of new histogram
    if dim == 2,
        [s,r] = strtok(h.type, '- ');
        hs.type = strtok(r, '- ');
    else
        hs.type = strtok(h.type, '- ');
    end
   
end
        