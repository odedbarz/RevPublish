function h = cat(dim, ha, hb, varargin)
% HISTOGRAM/CAT - Concatenate two or more histograms.
% H = CAT(DIM, HA, HB) concatenates the histograms HA and HB along dimension DIM.
% CAT works with arbitrary number of histograms and with neurograms.

% use recursion to handle case of more than 2 inputs
if nargin > 3,
    h = cat(dim, cat(dim, ha, hb), varargin{1}, varargin{2:end});
    return
end

dmax = ndims(ha.data);
if dim > dmax,
    error('Requested dimension out of range');
end

siza = size(ha.data);
sizb = size(hb.data);
siza(dim)=[];
sizb(dim)=[];
if length(siza) ~= length(sizb) | any(siza ~= sizb),
    error('Data sizes do not match');
end

if ~strcmp(ha.type, hb.type),
    error('Types do not match');
end
if any(ha.binwidth ~= hb.binwidth),
   error('Bin widths do not match');
end

nvar = length(ha.value);
if length(hb.value) ~= nvar,
    error('Stimulus variables do not match');
end
for k = 1:nvar,
    if k ~= dim & (length(ha.value{k}) ~= length(hb.value{k}) | ...
                   any(ha.value{k} ~= hb.value{k})),
        error(['Inconsistent ' ha.varname{k} ' vectors']);
    end
end

if dim > nvar & ~(dim == dmax-1 & ~isempty(strfind(ha.type, 'Raster'))),
    bwdim = dim - dmax + 2;
    a_end = ha.offset(bwdim) + size(ha.data,dim)*ha.binwidth(bwdim);
    if abs(hb.offset(bwdim) - a_end) > eps,
        error('Offsets do not match');
    end
end
        
% initialize result
h = ha;

% concatenate data
h.data = cat(dim, full(ha.data), full(hb.data));

% concatenate counts and stimulus variables
if dim <= nvar,
    h.synccount = cat(dim, ha.synccount, hb.synccount);
    h.spikecount = cat(dim, ha.spikecount, hb.spikecount);
    h.periodcount = cat(dim, ha.periodcount, hb.periodcount);
    h.value{dim} = [ha.value{dim}; hb.value{dim}];
elseif dim == dmax-1 & ~isempty(strfind(h.type, 'Raster')),
    h.synccount = ha.synccount + hb.synccount;  
    h.spikecount = ha.spikecount + hb.spikecount;
    h.periodcount = ha.periodcount + hb.periodcount;
end



    
