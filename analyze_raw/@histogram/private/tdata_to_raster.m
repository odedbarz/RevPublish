function data = tdata_to_raster(h, trialno, tdata)
% TDATA_TO_RASTER- Put time/trial data into 2D raster histogram
% tdata_to_raster(H, TRIALNO, TDATA) pust the trial # TRIALNO
% and time data TDATA into the bins of the 2D histogram H
% Spike and sync counts are not updated
%

% find column indices
colidx = 1+floor((tdata-h.offset(2))/h.binwidth(2));

h.offset(1) = 1;
h.binwidth(1) = 1;

% restrict indices to matrix dimensions
dmax = ndims(h.data);
ok = find(trialno >= 1 & trialno <= size(h.data,dmax-1) & ...
          colidx >= 1 & colidx <= size(h.data,dmax));

% use the Jay Delosge trick
data = sparse(trialno(ok),colidx(ok),1, size(h.data,dmax-1),size(h.data,dmax));
if nzmax(data) > prod(size(data))/4, 
    data = full(data); 
end

