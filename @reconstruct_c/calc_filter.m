function [G, pars] = calc_filter(obj, X_train, r_train, varargin)
%            
% function G = calc_filter(obj, X_train, r_train, [n_splits], [G], [lags], [bias])
% 


%% Set the input data
p = inputParser;

% Required
addRequired(p, 'X_train', @isnumeric);            
addRequired(p, 'r_train', @isnumeric);            

% Optional
addOptional(p, 'lags_ms', obj.lags_ms, @isnumeric);         % (ms)
addOptional(p, 'f', obj.f, @isnumeric);                     % (Hz) frequency axis     
addOptional(p, 'algo_type', obj.algo_type, @isstr); 
addOptional(p, 'iscausal', obj.iscausal, @(x) isnumeric(x) || islogical(x)); 
addOptional(p, 'bias', obj.bias, @isstr); 

addOptional(p, 'verbose', 0, @isnumeric);           
addOptional(p, 'fignum', [], @isnumeric);           

parse(p, X_train, r_train, varargin{:});
pars = p.Results;



%%
% Update the object, if needed
obj.lags_ms     = pars.lags_ms;     
obj.iscausal  	= pars.iscausal;
obj.algo_type 	= pars.algo_type;     
obj.bias        = pars.bias;

obj.X_train  	= X_train;
obj.r_train     = r_train;

if strcmpi('regression', obj.algo_type)
    % Apply ridge regression
    [obj.G, obj.Crr, obj.Crr_inv, obj.Crs, obj.eigv, obj.bias] = ...
        ridge_regression(X_train, r_train,...
        obj.lags,...
        obj.ridge.gamma, ...
        obj.ridge.inv_type,...
        obj.bias,...
        obj.ridge.xcorr_type,...
        obj.gpu_flag );    

    
elseif strcmpi('ASD', obj.algo_type)
    % BIAS
    if 1 == obj.bias.remove
        obj.bias.train.Sft  = mean(X_train, 2);    % spectrogram's bias 
        obj.bias.train.resp = mean(r_train);       % PSTH's bias
    else
        obj.bias.train.Sft  = 0.0;                 % spectrogram's bias 
        obj.bias.train.resp = 0.0;                 % PSTH's bias
    end

    % Subtract the bias before the reconstruction
    X_train = X_train - obj.bias.train.Sft;  
    r_train = r_train - obj.bias.train.resp;   

    
    % =======================================
    R       = response_mtx(r_train, obj.lags);
    x       = R';
    nks     = [obj.n_neurons*obj.n_lags, 1];  % number of filter pixels along [cols, rows]
    minlens = [2; 2];  % minimum length scale along each dimension

    obj.G = nan(obj.n_neurons*obj.n_lags, obj.n_bands);
    
    % Loop over all frequency bands seperately each filter
    Kasd = nan(obj.n_neurons*obj.n_lags, obj.n_bands);
    parfor fr = 1:obj.n_bands
        y = X_train(fr,:)';
        [kasd, ~] = fastASD(x, y, nks, minlens);            
        %obj.G(:,fr) = kasd;
        Kasd(:,fr) = kasd;
    end
    obj.G = Kasd;
    

elseif strcmpi('svd', obj.algo_type)
    % BIAS
    if 1 == obj.bias.remove
        obj.bias.train.Sft  = mean(X_train, 2);     % spectrogram's bias 
        obj.bias.train.resp = mean(r_train);        % PSTH's bias
        %obj.bias.train.r_std = std(r_train);        % variance
    else
        obj.bias.train.Sft  = 0.0;                  % spectrogram's bias 
        obj.bias.train.resp = 0.0;                  % PSTH's bias
    end

    % Subtract the bias before the reconstruction
    X_train = X_train - obj.bias.train.Sft;  
    r_train = r_train - obj.bias.train.resp;   
    %r_train = r_train./(eps + obj.bias.train.r_std);
    
    
    % =======================================
    R = response_mtx(r_train, obj.lags);

    % Calculate the reconstruction filters
    [U, S, V] = svds(R, size(R,1));
    
    % Use inv(S); S is a diagonal matrix 
    Sinv = inv(S);
    
    % Make sure that Sinv is well defined
    if min(diag(Sinv)) <= 1e3*eps        
        warning( '%s\n\t%s\n\t%s (gamma = %g)', ...
            '[reconstruct_c/calc_filter.m]', ... 
            '---> inv(S) might be ill-conditioned!',...
            '---> RIDGE REGRESSION was used',...
            obj.ridge.gamma );
        Sinv = diag(diag(Sinv) + obj.ridge.gamma);
        warning('--> [reconstruct_c/calc_filter.m]: inv(S) might be ill-conditioned!');
    end
    
    % R^{-T} = U*Sinv*V'
    obj.G = (U*Sinv*V')*X_train';
    
    
    
else
    error('--> [strf_c/calc_strf.m]: unrecognized ALGO_TYPE (%s) parameter!!!', pars.algo_type);
    
end















