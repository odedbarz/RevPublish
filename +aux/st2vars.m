function out = st2vars( st )
% 
% function varargout = st2vars( st )
% 
% Structure to variables
% 

if ~isstruct(st), return; end

fields = fieldnames(st);
% N = min(nargout, length(fields));
N = length(fields);
assert(nargout<=N, '--> Error in [st2vars.m]: There are too many output variables!!!')

out = cell(1,nargout);

for nn = 1:N
    out{nn} = st.(fields{nn});
end

