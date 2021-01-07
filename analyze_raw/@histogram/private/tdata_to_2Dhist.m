function data = tdata_to_2Dhist(h, rowdata, coldata)
% TDATA_TO_2DHIST- Put 2 sets of time data into 2D histogram
% TDATA_TO_2DHIST(H, ROWDATA, COLDATA) puts the 
% time data ROWDATA and COLDATA into the bins of the 2D histogram H.
% Spike and sync counts are not updated
%

% find row & column indices
rowidx = 1+floor((rowdata-h.offset(1))/h.binwidth(1));
colidx = 1+floor((coldata-h.offset(2))/h.binwidth(2));

% restrict indices to matrix dimensions
dmax = ndims(h.data);
ok = find(rowidx >= 1 & rowidx <= size(h.data,dmax-1) & ...
          colidx >= 1 & colidx <= size(h.data,dmax));

% use the Jay Delosge trick
data = sparse(rowidx(ok),colidx(ok),1, size(h.data,dmax-1),size(h.data,dmax));
if nzmax(data) > prod(size(data))/4, 
    data = full(data); 
end


