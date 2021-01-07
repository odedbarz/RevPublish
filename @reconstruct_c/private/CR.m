function Crs = CR(response, Sft, lags)
%
%   function Crr = AC(response, lags)
%
% Description:
% Implements,
%     Crs = R * Sft';
% But without calculating R.


% Calculate the autocorrelation matrix using for loops
n_units = size(response,2);
n_lags  = length(lags);
Crs     = nan(n_units*n_lags, size(Sft,1));  % initiate a square matrix

for n = 1:n_units
    Xn = lagged_mtx(response(:,n)', lags);
    idx_m = (n-1)*n_lags+(1:n_lags);
    Crs(idx_m,:) = Xn * Sft';
end
    
    
