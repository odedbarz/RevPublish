function y = MD(x, dim)
%
%   function y = MD(x, [dim])
%
% Modulation depth function (coefficient of variation).
%
% Description:
% The sqrt(2) coefficient of variation of the signal x. 

if 2 > nargin || isempty(dim)
    [~, dim] = max(size( x ));
end

% Modulation depth
y = sqrt(2) * std(x,[],dim)./mean(x,dim);
