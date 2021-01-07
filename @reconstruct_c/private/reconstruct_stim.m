function [Sest, pars] = reconstruct_stim(Sft, resp_train, resp_test,...
                            lags, inv_type, bias_choice, fignum)
%
%   function [Sest, pars] = reconstruct_stim(Sft, response_train, ...
%       response_test, lags, inv_type, bias_choice)
%
% Input:
%   Sft       : (bands x time) spectrogram spectrum (without the phase).
%   resp_train: (time x neurons) neural response for the training (building) 
%                of the reconstruction functions.
%   resp_test : (time x neurons) neural response for the testing of the
%                reconstruction functions.
%   lags      : (1 x K) the lags to calculate for G.
%   inv_type  : select the inversion type to use. 
%               inv_type == 0: (1x1, default) use MATLAB inversion operator /; 
%               inv_type == 1: (1x2) use "big" (threshold == inv_type(2)) 
%                              singular eigenvalues; 
%               inv_type == 2: singular eigenvalues with ridge regression 
%                              (ridge regression threshold == inv_type(2)); 
%
%   bias_choice: (1x2, cell) 
%           * {1}: this is a 1 or 0 (true or false) field; 
%           * {2}: this second entry, if given, is the bias to subtract from the
%                  estimated spectrogram. If it's not given (i.e., length(remove_bias) == 1) 
%                  then the bias will be computed from the training spectrogram.
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

if 3 > nargin
    resp_test = [];
end

% Choose the inversion type for the auto-correlation (AC) matrix
if 5 > nargin || isempty(inv_type)
    inv_type = 0;   % (default) use MATLAB's inversion operator /.
end

% Remove BIAS from spectrogram & responses (PSTHs) befor calculating the
% reconstruction filters
if 6 > nargin || isempty(bias_choice)
    bias.remove    = 1;             % (default) 
    bias.test.Sft  = [];            % (default) no bias to remove
elseif 6 <= nargin && 1 == length(bias_choice)
    bias.remove    = bias_choice(1);       % (default) 
    bias.test.Sft  = [];            % (default) no bias to remove
elseif 6 <= nargin && iscell(bias_choice) && 2 == length(bias_choice)
    bias.remove    = bias_choice{1};       % (default) 
    bias.test.Sft  = bias_choice{2};       % (default) no bias to remove
else
    error('--> [stim_reconstruction.m]: please choose something form the list above!');
end

% DEBUG-MODE
if 7 > nargin 
    fignum = [];   % (default) 
end

% Initialize the output variables
Sest = [];
G    = [];
eigv = [];


% # of lags for the reconstruction filters
n_lags = length(lags);
if isempty(n_lags), return; end



%% BIAS -- TRAINING -- REMOVE
if bias.remove
%     train.bias.Sft  = mean(train.Sft(:));        % spectrogram's bias 
    bias.train.Sft  = mean(Sft,2);              % spectrogram's bias 
    bias.train.resp = mean(resp_train);         % PSTH's bias
    
    % Subtract the bias
    Sft  = Sft - bias.train.Sft;  
    resp_train = resp_train - bias.train.resp;    
end



%% Build the auto-correlation (AC) of neuronl response
R = response_mtx(resp_train, lags);
Crr = R*R';



%% Build the cross-correlation stimulus & neural-response
% The i'th column holds the crosscorrelation inner product for the gi'th
%   reconstructed filter.
Crs = R * Sft';



%% Compute the inverse & RECONSTRUCTION filters
switch inv_type(1)
    case 0
        G = Crr \ Crs;
        
    case 1
        % Singular value threshold for the SVD decomposition
        singular_value_threshold = 1.0 - inv_type(2);
        
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
        singular_value_threshold = inv_type(2);
        
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



%% *** ESTIMATION ***
if ~isempty(resp_test)
    
    % BIAS -- TESTING -- REMOVE
    if bias.remove
        % Calc. the bias of the response
        bias.test.resp = mean(resp_test);	% PSTH's bias

        % Subtract the bias
        resp_test = resp_test - bias.test.resp;    
    end
    
    
    Rtest = response_mtx(resp_test, lags);

    % Calc. the estimated spectrogram
    Sest = G' * Rtest;
    
    % BIAS -- TESTING -- REMOVE
    if bias.remove
        if ~isempty(bias.test.Sft)  % does the user gave a bias?
            % Usually, this is the bias of the testing spectrogram to
            % compare with
            Sest = Sest + bias.test.Sft;
        else
            Sest = Sest + bias.train.Sft;
        end
    end
    
    
end


if 2 <= nargout
    pars.G       = G;
    pars.Crr     = Crr;
    pars.Crr_inv = Crr_inv;
    pars.Crs     = Crs;
    pars.eigv    = eigv;
    pars.bias    = bias;
end




%% Plot #1: the reconstruction
if isempty(fignum), return; end

figure(fignum);
clf;

ax = subplot(2,1,1);
imagesc(Sft);
ylabel('Target S(f,t)');

ax(2) = subplot(2,1,2);
imagesc(Sest);
ylabel('$\hat{S}$(f,t)');

linkaxes(ax);







