function n = num_dim(type)
% NUM_DIM - Return expected number of dimensions for histogram type
% NUM_DIM(type) returns 1 or 2 depending on whether histogram is 1D or 2D
%
if strcmp(type, 'PESE'),
     n = 2;
elseif ~isempty(find(type == ' ' | type == '-')),
     n = 2;
else n = 1;
end
