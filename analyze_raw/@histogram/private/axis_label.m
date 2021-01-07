function s = axis_label(h, dim)
% AXIS_LABEL - Get axis label from histogram type
% AXIS_LABEL(h, DIM) returns a label for dimension DIM in H.
% If dim is missing, the default dimension is returned for 1D histograms.
% AXIS_LABEL(H, what), where what is 'rate', or 'count' returns a label
% for the bin counts.
%

if nargin < 2, dim = which_dim(h); end

if ischar(dim),

    if ~isempty(strfind(h.type, 'Synchrony')),
        s = 'Synchrony';
    elseif strcmp(dim, 'rate'),
        s = 'Discharge Rate (sp/s)';
    elseif ~isempty(strfind(h.type, 'correlation')),
        s = 'Number of Intervals';
    else
        s = 'Number of Spikes';
    end

else    % numeric dimension

    if dim <= length(h.value),  % neurogram variable

        if dim <= length(h.varunit) & ~isempty(h.varunit{dim}),
            s = [h.varname{dim} ' (' h.varunit{dim} ')'];
        else
            s = h.varname{dim};
        end

    else        % histogram variable

        if num_dim(h.type) == 2 & dim == ndims(h.data)-1,   % Y dimension comes second in type
            [tmp, htype] = strtok(h.type, ' -');
        else htype = h.type;
        end

        s = type2label(strtok(htype, ' -'));
    end

end

