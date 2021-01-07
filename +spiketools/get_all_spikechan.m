function all_spikechan = get_all_spikechan(S, syncchan)
%
%   function all_spikechan = get_all_spikechan(S, syncchan)
%
% Description:
% Get all spike channels (spikechan) from S. S can be an Impale structure
% or a structure array.
%

if 2>nargin, syncchan = 0; end

all_spikechan = cell(size(S));

for kk = 1:length(S)
    if iscell(S)
        Sk = S{kk}.ch;
    else
        Sk = S(kk).ch;
    end
    dummy = cellfun(@(C) C(C~=syncchan)', Sk, 'UniformOutput', false);
    all_spikechan{kk} = unique( [dummy{:}] );
end

all_spikechan = unique( cell2mat(all_spikechan) );
