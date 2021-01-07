function data = double(h, type)
% HISTOGRAM/DOUBLE - Convert histogram to double
% data = double(h) returns the bin counts in h.
% By default, the data are returned as a full matrix even if they are
% stored in the histogram as a sparse matrix.  Data = double(h, 'raw')
% returns the data exactly as they are stored in the histogram. 
%
if nargin > 1 & strcmp(type, 'raw'),
    data = h.data;
else
    data = full(h.data);
end
