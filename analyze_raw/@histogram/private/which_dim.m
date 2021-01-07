function dim = which_dim(h)
% WHICH_DIM - Find best dimension for histogram
% Normally, it is the second (column) dimension
%
dim = max(find(size(h.data)~=1));
if isempty(dim), dim = ndims(h.data); end