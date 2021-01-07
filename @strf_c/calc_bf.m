function [BF, amp, row, col] = calc_bf(obj)
%
% function [BF, amp, row, col] = calc_bf(obj)
%
% Description: 
% The best frequency (BF) of the STRF is its maximum point.
%

[n_bands, n_lags] = size(obj.strf);
assert(n_bands == length(obj.f));

% The BF is the max frequency (y-axis) in the STRF 
[amp, idx] = max( obj.strf(:) );
[row, col] = ind2sub([n_bands, n_lags], idx);
BF = obj.f(row);   % (Hz) 


