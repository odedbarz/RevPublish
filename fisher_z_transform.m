function z = fisher_z_transform(r)
%
%   function z = fisher_z_transform(r)
%
% Description:
%
%   z = atanh(r)
%
%   In statistics, the Fisher transformation (aka Fisher z-transformation) can 
%   be used to test hypotheses about the value of the population correlation 
%   coefficient ρ between variables X and Y.[1][2] This is because, when the 
%   transformation is applied to the sample correlation coefficient, the sampling 
%   distribution of the resulting variable is approximately normal, with a variance 
%   that is stable over different values of the underlying true correlation. 
%
%   Source: Wikipedia (https://en.wikipedia.org/wiki/Fisher_transformation)
%

z = atanh(r); 
