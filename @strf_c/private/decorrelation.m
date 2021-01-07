function [strf, strf_noncausal, Clags] = decorrelation(Cw, Rw, gamma, iscausal, verbose)
%
%   [strf, strf_noncausal, Clags] = decorrelation(Cw, Rw, [tol], [iscausal], [verbose])
%
%

%% Check the input
% Ridge regression tolerance, as a ratio of MAX_SINGULAR_VALUE
if 3 > nargin || isempty(gamma)
    gamma = 0.001;
end

if 4 > nargin || isempty(iscausal)
    iscausal = 1;
end

if 5 > nargin || isempty(verbose)
    verbose = 0;
end


%%
[n_bands, n_lags] = size(Rw);

% Inverse of the DFT Autocorrelation Matrix (AC) 
%n_lags = 2*n_win + 1;
assert( (n_lags - 1)/2 == fix((n_lags - 1)/2),...
    '--> ERROR at [decorrelation.m]: set the size of the temporal window to be of odd length!' )
n_win = (n_lags - 1)/2;

% The STRF's temporal lags
lags = -n_win:n_win;

% off-(lower)diagonal and diagonal indices
n_tril = (n_bands+1)*n_bands/2;     


% Construct the indices that assign the rows in Aavg to a matrix At in the right order
col2mtx_idx = zeros(n_bands, n_bands);
assert(size(Cw,1) == n_tril,...
    '--> ERROR @ (size(Cw,1) == n_tril): these two MUST be the same!!!');
col2mtx_idx( 1 == tril(ones(n_bands)) ) = 1:n_tril;
col2mtx_idx = col2mtx_idx + tril(col2mtx_idx,-1)';     % symmetric indices for a symmetric auto-correlation matrix

% Cell array to hold all U, S, & V matrices 
U = cell(1, n_win+1);
D = cell(1, n_win+1);
V = cell(1, n_win+1);

% Get the maximum singular value from all lag matrices. This value will
% be used to find the best ridge regression parameter (via jackknife).
max_singular_value = -1;

Clags = cell(1, n_win+1);

for k = 1:(n_win+1)     % n_win+1, +1 for the 'middle' lag
    % From column to an AC matrix
    Clags{k} = reshape(Cw(col2mtx_idx, k), n_bands, n_bands);    
    
    [U{k}, Dk, V{k}] = svd(Clags{k});
    D{k} = diag(Dk);
    
    % Get the maximum eigenvalue of all autocorrelation matrices
    max_singular_value = max([max_singular_value, D{k}(1)]);
end

if verbose
    fprintf('--> max_singular_value: %g \n', max_singular_value);
end



%% STRF; Ridge-regression
% Invert the the AC matrix, and for each lag, multiply by the stimulus-response
% cross-correlation vector. Then, sum it all up to get the STRF

% The STRF, in the frequency domain
strf_w = nan(n_bands, n_lags);

for k = 1:(n_win+1)     % n_win+1, +1 for the 'middle' lag
    % The inverse eigenvalue + ridge regression
    Dinv = 1./(gamma*max_singular_value + D{k});
    
    % Construct the STRF 
    strf_w(:,k) = (V{k} * diag(Dinv) * U{k}') * Rw(:,k);
    
end

% Force the STRF(w) to be conjugate symmetric along its temporal axis (that is, 
% make the STRF real in the time domain)
strf_w = [strf_w(:,1:(n_win+1)), conj( strf_w(:,(1+n_win):-1:2) )];

strf_noncausal = ifft(strf_w, [], 2, 'symmetric');

% Causal temporal window (following Theunissen @ STRFLAB)
if 1 == iscausal
    % Create a "causal window" to dispose the negative lags (STRFLAB)
    %win_causal = (atan(lags)+pi/2)/pi;  
    %strf_causal = strf_noncausal .* win_causal;
    %strf = strf_causal(:, lags>=0);
    
    % Don't use the "causal" window
    strf = strf_noncausal(:, lags>=0);  
    
else
    strf = strf_noncausal;
    
end


if 0 ~= nnz(sum(imag(ifft(strf_w, [], 2))))
    fprintf('--> [strf_w]: STRF(f,lags) is not conjugate symmetric along the temporal axis!!!\n');
end










