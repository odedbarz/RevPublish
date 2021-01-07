function r_est = predict(obj, X_test, varargin)
%
% function r_est = predict(X_test, varargin)
%
%


%% Set the inputs
p = inputParser;
addRequired(p, 'obj', @isobject);       % strf object
addRequired(p, 'X_test', @isnumeric);	% (n_bands x nt) spectrogram

% % (n_time x 1) split the spectrogram into different stimuli
% addOptional(p, 'chunk_idx', ones(1,size(X_test,2)), @isnumeric);       

% addOptional(p, 'strf', obj.strf, @isnumeric);     % (n_bands x n_lags) strf
addOptional(p, 'avg', [], @isnumeric);           
addOptional(p, 'normalize', [], @isnumeric);           
parse(p, obj, X_test, varargin{:});
pars      = p.Results;

normalize = pars.normalize;
avg       = pars.avg;

assert(~isempty(obj.strf), '--> [@strf_c\predict.m]: You need first to FIT a STRF!');


%% Estimate response using STRF
% CONVOLUTION: Option #1
%{
% n_smp  : # of samples 
% n_bands: # of frequency bands
[n_bands, n_smp] = size(Sft);

conv_Sft_by_strf = arrayfun(@(I) conv(Sft(I, :)', strf(I,:)), 1:n_bands, 'UniformOutput', 0);
rLin = [conv_Sft_by_strf{:}]';

% Remove extra samples, if needed. This makes sure that the linear response 
% PSTH has the same length as that of the spectrogram. Sometimes, due to
% round off errors, there is a discrepancy between the two of few samples.
rLin = sum(rLin(:, 1:n_smp), 1);

% Make the estimated response a column vector
rLin = rLin(:);
%}

% CONVOLUTION: Option #2
% Xconv = convmtx_spec(X_test, size(obj.strf,2));
Xconv = convmtx_spec(X_test, obj.lags);
r_est = Xconv' * obj.strf(:);

% Add back the means
if ~isempty(avg)
    r_est = r_est + avg;
end

if ~isempty(normalize)
    r_est = r_est/max(r_est) * normalize;
end

obj.r_est = r_est;



