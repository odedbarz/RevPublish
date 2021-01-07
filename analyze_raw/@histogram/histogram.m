function h = histogram(data)
% HISTOGRAM - Histogram class constructor
% H = HISTOGRAM(DATA) creates a histogram whose bin counts are given by DATA.
% If DATA has N > 2 dimensions, a neurogram with N-2 stimulus variables is created.
%
if nargin == 0 || isstruct(data),
    h.data = 0;
elseif isa(data, 'histogram'),
    h = data;
    return
 elseif ~isstruct(data),
    h.data = data;
end

h.type = 'none';
h.binwidth = [1 1];
h.offset = [0 0];
h.syncchan = 0;
h.spikechan = 1;
h.period = NaN;
h.order = 1;
h.stimfile = '';

% counts to be filled in
if ndims(h.data) > 2,   % neurogram
    csize = size(h.data);
    csize = csize(1:end-2);
    if length(csize) < 2, csize = [csize 1]; end
    h.spikecount = zeros(csize);
    h.synccount = zeros(csize);
    h.periodcount = zeros(csize);
else
    h.spikecount = 0;
    h.synccount = 0;
    h.periodcount = 0;
end

% gates
g.type = 'none';
g.delay = 0;
g.width = Inf;
h.gate = g;

% for neurograms
h.varname = {};
h.value = {};
h.varunit = {};
for dim = 1:ndims(h.data)-2,
    h.value{dim} = [0:size(h.data,dim)-1]';
    h.varname{dim} = ['Variable' num2str(dim)];
end

h.userdata = [];

if nargin>0 && isstruct(data),
   flds = fieldnames(data);
   for kf = 1:length(flds),
      if isfield(h, flds{kf}),
         h.(flds{kf}) = data.(flds{kf});
      end
   end
   if isfield(data, 'size'),
      h.data = zeros(data.size);
   end
end


h = class(h, 'histogram');

