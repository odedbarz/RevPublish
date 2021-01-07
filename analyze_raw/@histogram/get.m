function val = get(h, prop_name, dim)
% HISTOGRAM/GET - Get values of histogram properties
% Valid properties are 'Type', 'Size', 'BinWidth', 'Offset', 'BinTimes',
%   'SyncCount', 'SpikeCount', 'PeriodCount', 'SyncChan', 'SpikeChan',
%   'GateType', 'GateDelay', 'GateWidth', and 'UserData'
% Additional properties for neurograms are 'VarName', 'VarUnit', and any
% named stimulus variable
%
if nargin < 3, dim = which_dim(h); end

switch prop_name
   case 'Type'
      val = h.type;
   case 'NumDim'
      val = num_dim(h.type);
   case 'Size'
      val = size(h.data);

   case 'BinWidths'
      val = h.binwidth;
   case 'BinWidth'
      val = h.binwidth(dim);
   case 'Offsets'
      val = h.offset;
   case 'Offset'
      val = h.offset(dim);
   case 'BinTimes'
      bwdim = dim - ndims(h.data) + 2;
      val = h.offset(bwdim) + [0:size(h.data,dim)-1]*h.binwidth(bwdim);

   case 'SyncChan'
      val = h.syncchan;
   case 'SpikeChan'
      val = h.spikechan;

   case 'Period'
      if isempty(strfind(h.type, 'Period')),
         val = NaN;
      else val = h.period;
      end
   case 'Order'
      val = h.order;
   case 'StimFile'
      if strcmp(h.type, 'Revcor') | strcmp(h.type, 'PESE'),
         val = h.stimfile;
      else val = '';
      end

   case 'GateType'
      val = h.gate.type;
   case 'GateDelay'
      val = h.gate.delay;
   case 'GateWidth'
      val = h.gate.width;

   case 'SyncCount'
      val = h.synccount;
   case 'SpikeCount'
      val = h.spikecount;
   case 'PeriodCount'
      val = h.periodcount;

   case 'VarName'
      val = h.varname{dim};
   case 'VarUnit'
      val = h.varunit{dim};
   case 'NumVar'
      val = length(h.value);

   case 'UserData'
      val = h.userdata;

   otherwise
      dim = strmatch(prop_name, h.varname, 'exact');
      if ~isempty(dim)
         val = h.value{dim};
      else
         error([prop_name ' is not a valid histogram property']);
      end
end
