function [X_est, R_test] = predict(obj, r, varargin)
%    
% function [X_est, R_test] = predict(obj, r, [G], [lags], [bias])
% 
% Input:
%   obj : reconstruction object
%   r   : response matrix to use for the test estimation
% 
%

%% Set the input data
p = inputParser;
addRequired(p, 'resp_test', @isnumeric);            % (n_bands x n_time) spectrogram power

addOptional(p, 'G', obj.G, @isnumeric);             % (ms)
addOptional(p, 'lags', obj.lags, @isnumeric);           
addOptional(p, 'bias', obj.bias, @isstruct);           
addOptional(p, 'fignum', [], @isnumeric);           

parse(p, r, varargin{:});
pars = p.Results;

G       = p.Results.G;       
lags    = p.Results.lags;       
bias    = p.Results.bias;       

assert( size(r,2) == obj.n_neurons,...
    sprintf('--> [reconstruct_c/predic.m]:\n\t%s\n\t%s', ...
            'the number of columns in the response matrix is wrong!',...
    sprintf('The response matrix r should have %d columns/responses\n', obj.n_neurons)) );




%%

% BIAS -- TESTING -- REMOVE
if bias.remove
    % Calc. the bias of the response
    bias.test.resp = mean(r);	% PSTH's bias
    
    % Subtract the bias
    r = r - bias.test.resp;  
    
    %{
    if strcmpi('svd', obj.algo_type)
        bias.test.r_std= std(r);
        r = r./(eps + bias.test.r_std);
    end
    %}
    
    % ZCA
    %r = (r - bias.test.resp)./(eps + bias.test.r_std);
    
end


% R_test is a convolution matrix such that R_test' * g, and g is a column
% in the filter matrix G
R_test = response_mtx(r, lags);

% Calc. the estimated spectrogram
X_est = G' * R_test;

% BIAS -- TESTING
if bias.remove
    % Add back the REMOVED bias
    X_est = X_est + bias.train.Sft;
end

% Save into the object for later use
obj.X_est  = X_est;
obj.R_test = R_test;


%% Plot the reconstructed spectrogram
if isempty(pars.fignum) || pars.fignum<=0, return; end

figure(pars.fignum);
clf;

imagesc(X_est);
set(gca, 'YDir', 'normal');
title('$\hat{S}$(f,t)');



