function str2 = mName2latex(str)
%
%   function str2 = mName2latex(str)
%
% Description:
%   Replaces all underscores ('_') into '\_' for the Latex interpretor

if iscell(str)
    non_empty = cellfun(@(X) ~isempty(X), str);
    str2 = cellfun(@(S) replace(S, '_', '\_'), str(non_empty), 'UniformOutput', false);
    str2 = cellfun(@(S) replace(S, '%', '\%'), str2, 'UniformOutput', false);
    str2 = cellfun(@(S) replace(S, '#', '\#'), str2, 'UniformOutput', false);
else
    str2 = replace(str, '_', '\_');
    str2 = replace(str2, '%', '\%');
end

