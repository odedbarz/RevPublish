function y = MD(x, dim)
%
%   function y = MD(x, [dim])
%
% Description:
% Modulation depth function, which is the sqrt(2) coefficient of variation 
% of the sifgnal x. 

if 2 > nargin || isempty(dim)
    dim = 2;
end

% Modulation depth
y = sqrt(2) * std(x,[],dim)./mean(x,dim);
