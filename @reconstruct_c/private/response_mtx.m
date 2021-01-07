function R = response_mtx(response, lags, algo_type)
%
% function R = response_mtSft(response, lags, [algo_type])
%
% Input:
%   response: (n_time x n_neurons) spectrogram.
%   lags    : (1 x n) a vector of temporal lags, starting at zero.
%   [algo_type]:
%       algo_type == 1: use convmtx.
%       algo_type == 2: use for loop (code from Nima Mesgarani lab).
%
% Description:
% Create a response\convolution matriSft.
%
% * n_time   : # of time samples
% * n_neurons: # of neurons in the training
%
% (OBZ) Notes:
% See github, from Nima Mesgarani's lab:
%   https://github.com/odedbarz/Naplib/blob/master/Matlab/Reconstrcution/LagGeneratorNew.m
%
% The other option is to use a CONVMTX instead of this for loop. The results 
% are identical up to a transform on the output matrices (in my function, 
% I flaten the matrix so you need to reshape it first). 
% 
% For example:
% >> A = LagGeneratorNew((1:5)', -3:3)' 
% >> B = response_mtx((1:5)', -3:3)
% >> isequal(A, B)  % ==> 1 
%
 

[n_time, n_neurons] = size(response);

% # of lags for the reconstruction filters
n_lags = length(lags);
if isempty(n_lags), return; end

if 2 >= nargin
    algo_type = 1;
end

%%
R = nan(n_lags*n_neurons, n_time, 'like', response);

switch algo_type
    case 1
        % Handle positive delays, if needed (lags(1) > 0)
        if 0 < lags(1)
            response = [zeros(lags(1), size(response,2)); response(1:(end-lags(1)),:)];
            lags = lags - lags(1);
        end

        for k = 1:n_neurons
            Xk = convmtx(response(:,k)', n_lags);
            Xk = Xk(:,(1-lags(1)):(end-lags(end)));
            R((k-1)*n_lags+(1:n_lags), :) = Xk;
        end
        
    case 2
        for k = 1:n_neurons
            Xk = LagGeneratorNew(response(:,k), lags)';
            R((k-1)*n_lags+(1:n_lags), :) = Xk;
        end

    otherwise
        error('--> [response_mtx.m]: Unrecognized ALGO_TYPE parameter!');
end




