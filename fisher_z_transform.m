function y = fisher_z_transform(x, trans_type)
%
%   function y = fisher_z_transform(x, trans_type)
%
%
%   trans_type: (Optional) if 'direct' (default) then perform the Fisher transform. 
%             If 'inv', transform the INVERSE transform. 
%
% Description:
%
%   z = atanh(r)
%
%   In statistics, the Fisher transformation (aka Fisher z-transformation) can 
%   be used to test hypotheses about the value of the population correlation 
%   coefficient Ï? between variables X and Y.[1][2] This is because, when the 
%   transformation is applied to the sample correlation coefficient, the sampling 
%   distribution of the resulting variable is approximately normal, with a variance 
%   that is stable over different values of the underlying true correlation. 
%
%   Source: Wikipedia (https://en.wikipedia.org/wiki/Fisher_transformation)
%

if 2 > nargin
    trans_type = 'direct';
end

if strcmpi(trans_type, 'direct')
    % z = atanh(r)
    y = atanh(x);
elseif strcmpi(trans_type, 'inv')
    % r = atanh(z)   
    y = tanh(x);
end
