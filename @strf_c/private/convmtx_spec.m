function X = convmtx_spec(Sft, lags)
%
% function X = convmtx_spec(Sft, lags)
%
% Input:
%   Sft : (n_time x n_neurons) spectrogram.
%   lags: (1 x n) a vector of temporal lags, starting at zero.
%
% Description:
% Create a response\convolution matrix.
%

[n_bands, n_time] = size(Sft);
n_lags = length(lags);

%{
% Didn't check this block
if 0 < lags(1)
    %response = [zeros(lags(1), size(response,2)); response(1:(end-lags(1)),:)];
    Sft = [zeros(size(Sft,1), lags(1)), Sft(:, 1:(end-lags(1)))];
    lags = lags - lags(1);
end
%}

X = zeros(n_lags*n_bands, n_time);
n0 = 1-lags(1);
for k = 1:n_lags
    %X((k-1)*n_bands+(1:n_bands), k:end) = Sft(:, 1:end-k+1);
    
    X((k-1)*n_bands+(1:n_bands), k:(end-n0+1)) = Sft(:, n0:end-k+1);
end


% 'ERROR!!!'
% '---> [STRF_C/convmtx_spec]: !!! check the X indices: (k-1)*n_bands+(1:n_bands) !!!'
