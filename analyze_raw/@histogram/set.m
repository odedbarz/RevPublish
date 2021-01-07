function h = set(h, varargin)
% HISTOGRAM/SET - Set values of histogram properties
% Valid properties: 'Type', 'Size', 'BinWidth(s)', 'Offset(s)', 'Period', 'Order'
% 'SyncChan', 'SpikeChan', 'GateType', 'GateDelay', 'GateWidth', 'StimFile', and 'UserData'
% Additional properties for neurograms are 'VarName(s)', 'VarUnit(s)', and any
% previously named stimulus variable
%
% Type histogram_types' to see all histogram types.
%
prop_list = varargin;

while length(prop_list) >= 2,
    
    prop_name = prop_list{1};
    val = prop_list{2};
    if length(prop_list) >= 3 && isnumeric(prop_list{3}),
            dim = prop_list{3};
            prop_list = prop_list(4:end);
    else
            dim = which_dim(h);
            prop_list = prop_list(3:end);
    end    
    
    switch prop_name
    case 'Type'
        check_type(val);
        h.type = val;
    case 'Size'
        if length(val) < num_dim(h.type),
             error('Size must be 1 X 2 vector');
        end
        h = resize(h, val);
        
    case 'BinWidths'
        if length(val) > 2 || length(val) < num_dim(h.type),
            error('Bin width must be 1 X 2 vector');
        end
        if any(val <= 0), error('Bin widths must be > 0'); end
        if length(val) == 2,
             h.binwidth = val;
        else h.binwidth(which_dim(h)) = val;
        end
    case 'BinWidth'
        if val <= 0, error('Bin widths must be > 0'); end
        h.binwidth(dim) = val;
        
    case 'Offsets'
        if length(val) > 2 || length(val) < num_dim(h.type),
            error('Offset must be 1 X 2 vector');
        end
        if length(val) == 2,
             h.offset = val;
        else h.offset(which_dim(h)) = val;
        end
    case 'Offset'
        h.offset(dim) = val;
        
    case 'SyncChan'
        h.syncchan = val;
    case 'SpikeChan'
        h.spikechan = val;
        
    case 'Period'
        if val <= 0, error('Period must be > 0'); end
        h.period = val;
    case 'Order'    % order of interval histogram
        if val <= 0, error('Order must be > 0'); end
        h.order = val;       
    case 'StimFile'
        h.stimfile = val;
        
    case 'GateType'
        switch val
        case 'none'
            h.gate.type = val;
            h.gate.delay = 0;
            h.gate.width = Inf;
        case {'PST', 'time', 'trial' }
            h.gate.type = val;
        otherwise
            error(['Illegal gate type: ' val]);
        end
    case 'GateDelay'
        h.gate.delay = val;
    case 'GateWidth'
        h.gate.width = val;
        
    case 'VarNames'
        if ~iscell(val),
            error('List of variable names must be cell array of strings')
        end
        h.varname = val;
        if length(h.value) > length(val),
           h.value = h.value(1:length(val));
        end
    case 'VarName'
        h.varname{dim} = val;
    case 'VarUnits'
        h.varunit = val;
    case 'VarUnit'
        h.varunit{dim} = val;
        
    case 'UserData'
        h.userdata = val;
        
    otherwise
        dim = strmatch(prop_name, h.varname, 'exact');
        if ~isempty(dim)
            h = set_stim_var(h, val(:), dim);
        else
            error([prop_name ' is not a valid histogram property']);
        end
    end %case
end % while 

%----------------------------------------------------------------
function h = set_stim_var(h, val, dim)
% Set stimulus variable in neurogram
%
if any(h.data(:) ~= 0),     % empty histogram: resize
    error('Neurogram must be empty to set variable values');
end

h.value{dim} = val;

% resize histogram to match variable length
old_size = size(h.data);
if dim <= length(old_size) - num_dim(h.type),   % replacing existing stimulus variable
    if numel(val) == old_size(dim), return; end     % size already correct
    new_size = [old_size(1:dim-1) numel(val) old_size(dim+1:end)];
else  % new stimulus variable
    new_size = [old_size(1:dim-1) numel(val) old_size(end-num_dim(h.type)+1:end)];
end

h = resize(h, new_size);


    
        
