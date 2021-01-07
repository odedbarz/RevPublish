function Xv = lagged_mtx(v, lags)
%
%   function X = lagged_mtx(v, lags)
%

if iscolumn(v)
    v = v';
elseif ~iscolumn(v) && ~isrow(v)
    error('--> [lagged_mtx.m] v must be a vector!');    
end

n_lags = length(lags);

Xv = convmtx(v, n_lags);
Xv = Xv(:,(1-lags(1)):(end-lags(end)));

