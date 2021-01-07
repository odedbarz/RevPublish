function check_type(type)
% CHECK_TYPE - Check histogram type
%
if ~isstr(type)
    error('Histogram type must be a string')
end
types = histogram_types;
if isempty(strmatch(type, types ,'exact')),
    error([type ' is an illegal histogram type'])
end
