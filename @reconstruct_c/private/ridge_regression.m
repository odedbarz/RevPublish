function [G, Crr, Crr_inv, Crs, eigv, bias] = ridge_regression(Sft, response, lags,...
    gamma, inv_type, bias, ...
    xcorr_type, gpu_flag)

%
%   [G, Crr, Crr_inv, Crs, eigv, bias] = ...
%       fit_ridge(Sft, response, lags, inv_type, bias_choice,...
%       [algo_type], [gpu_flag])
%
% Input:
%   Sft       : (bands x time) spectrogram spectrum (without the phase).
%   resp_train: (time x neurons) neural response for the training (building) 
%                of the reconstruction functions.
%   lags      : (1 x K) the lags to calculate for G.
%   inv_type  : select the inversion type to use. 
%               inv_type == 0: (1x1, default) use MATLAB inversion operator /; 
%               inv_type == 1: (1x2) use "big" (threshold == gamma) 
%                              singular eigenvalues; 
%               inv_type == 2: singular eigenvalues with ridge regression 
%                              (ridge regression threshold == gamma); 
%
%   bias_choice: (1x1) this is a 1 or 0 (true or false) field; 
%   fignum    : (1x1) figure number. If specified, the function will plot 
%                some of the results for debug purpose. Otherwise, set
%                as empty.
% Output:
%   Sest: (bands x time) estimated response
%   pars: (struct)
%       G   : ([lags*n_neurons] x bands) the reconstruction filter, one for each frequency
%             band.
%       (eigv): (struct) singular value data.
%       (bias): (struct) removed biased values
%       Crr, Crs, Crr_inv
%
% Reference:
%   Mesgarani, N., David, S. V., Fritz, J. B., & Shamma, S. A. (2009). Influence 
%   of context and behavior on stimulus reconstruction from neural activity in 
%   primary auditory cortex. Journal of Neurophysiology, 102(6), 3329–3339. 
%   https://doi.org/10.1152/jn.91128.2008
%
%   See also Bahar Khalighinejad (baharkh22) Github page at
%   https://github.com/Naplib/Naplib/tree/master/Matlab/Reconstrcution
%

if ~exist('xcorr_type', 'var')
    xcorr_type = 1;
end

if ~exist('gpu_flag', 'var')
    gpu_flag = 0;
elseif (1 == gpu_flag) && (1 == xcorr_type || 3 == xcorr_type)
    gpu_flag = 0;
    warning('\n\t%s\n\t%s\n',...
        '--> [fit_ridge.m]: GPU_FLAG is ON but you don''t have a GPU on this platform!',...
        '--> [fit_ridge.m]: turning off the GPU_FLAG!');
end

% Choose the inversion type for the auto-correlation (AC) matrix
if 4 > nargin || isempty(inv_type)
    inv_type = 2;   % (default) use MATLAB's inversion operator /.
end

% Initialize the output variables
G    = [];
eigv = [];

% # of lags for the reconstruction filters
n_lags = length(lags);
if isempty(n_lags), return; end



%% BIAS -- TRAINING -- REMOVE
if 1 == bias.remove
    bias.train.Sft  = mean(Sft, 2);         % spectrogram's bias 
    bias.train.resp = mean(response);       % PSTH's bias
    %bias.train.r_std = std(response);
else
    bias.train.Sft  = 0.0;                  % spectrogram's bias 
    bias.train.resp = 0.0;                  % PSTH's bias
end

% Subtract the bias before the reconstruction
Sft = Sft - bias.train.Sft;  

response = response - bias.train.resp;     
%response = response./bias.train.r_std;


% Use the GPU
if gpu_flag
    response = gpuArray(response);
end
    



%% Calculate correlations
if (1 == xcorr_type || 2 == xcorr_type)
    % Build the auto-correlation (AC) of neuronl response
    R   = response_mtx(response, lags, xcorr_type);
    Crr = R*R';
    
    % Build the cross-correlation stimulus & neural-response
    % The i'th column holds the crosscorrelation inner product for the gi'th
    %   reconstructed filter.
    Crs = R * Sft';   
    
elseif 3 == xcorr_type
    Crr = AC(response, lags);           % autocorrelation matrix
    Crs = CR(response, Sft, lags);      % cross-correlation matrix 
    
else
    error('--> ridge_regression.m: unrecognized ALGO_TYPE!');
    
end



%% Compute the inverse & RECONSTRUCTION filters
switch inv_type(1)
    case 0
        G = Crr \ Crs;
        Crr_inv = inv(Crr);
        
        
    case 1
        % Singular value threshold for the SVD decomposition
        singular_value_threshold = 1.0 - gamma;
        
        [U,S,V] = svd(Crr, 'econ');

        s = diag(S);

        % Finds the indices to invert for the SVD decomposition
        idx_s = cumsum(s)/sum(s) <= singular_value_threshold;
        n_s   = nnz(idx_s);     % # of non-zero singular values 
        s_inv = 1./s(idx_s);
        
        % Save for later use
        eigv.s = s;
        eigv.idx_s = idx_s;
        
        Crr_inv = V(:,1:n_s) * diag(s_inv) * U(:,1:n_s)';

        % The reconstruction filter(s)
        G = Crr_inv * Crs;
        

    case 2  % Ridge regression
        % Singular value threshold for the SVD decomposition
        singular_value_threshold = gamma;
        
        [U,S,V] = svd(Crr, 'econ');

        s = diag(S);
        s_inv = 1./(s + max(s)*singular_value_threshold);
        
        % Save for later use
        eigv.s = s;
        eigv.ridge = singular_value_threshold;
        
        Crr_inv = V * diag(s_inv) * U';

        % The reconstruction filter(s)
        G = Crr_inv * Crs;
        

    otherwise
        error('--> Error at [stim_reconstruction.m]: unrecognized INV_TYPE (inv_type = %g)!\n', inv_type);
end


%%
% Back to CPU
if gpu_flag
    G = gather(G);
    Crr = gather(Crr);
    Crr_inv = gather(Crr_inv);
    Crs = gather(Crs);
    eigv.s = gather(eigv.s);
    bias = gather(bias);
end










