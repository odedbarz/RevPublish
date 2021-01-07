function Crr = AC(response, lags)
%
%   function Crr = AC(response, lags)
%
% Description:
% Implements,
%     R   = response_mtx(response, lags, algo_type);
%     Crr = R*R';
% But without calculating R.

% Calculate the autocorrelation matrix using for loops
n_units = size(response,2);
n_lags  = length(lags);
Crr     = nan(n_units*n_lags);  % initiate a square matrix

for n = 1:n_units
    Xn = lagged_mtx(response(:,n)', lags);
    idx_n = (n-1)*n_lags+(1:n_lags);

    for m = 1:n_units
        Xm = lagged_mtx(response(:,m)', lags);
        idx_m = (m-1)*n_lags+(1:n_lags);
        Crr(idx_n, idx_m) = Xn * Xm';
    end
    
end
    
