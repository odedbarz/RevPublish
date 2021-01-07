function [Cx, bias] = xcorr_avg(X, Y, n_win, xcorr_type, remove_bias, gpu_flag)
%
% function [Cx, bias] = xcorr_avg(X, Y, n_win, xcorr_type, remove_bias, obj.gpu_flag)
%
% 
% Description:
% Returns the averaged autocorrelation matrix, <X(t),Y*(t)>, in the time
% domain. This matrix will later on be transformed into the frequency
% domain (<X(w),Y*(w)>; see Theunissen et al., 2000).
%

if isempty(Y)
    % Perform autocorrelation
    Y = X;
end

assert( max(size(X)) == max(size(Y)),...
    '--> ERROR at [xcorr_avg.m]: both X and Y MUST share same MAX dimension!' );
if isvector(Y) 
    Y = Y(:)';  % set Y to be a row vector     
end

if 4 > nargin || isempty(xcorr_type)
    xcorr_type = 1;
end

if 5 > nargin || isempty(remove_bias)
    remove_bias = false;     % remove the bias (unbias) from the cross-correlation
end

if 6 > nargin || isempty(gpu_flag)
    gpu_flag = false;     
end


%% ### GPU ###
if ~gpu_flag
    variable_type = 'double';
else
    variable_type = 'gpuArray';
end


%% Set the input variables
% Correlation parameters
maxlag = n_win;      	% (samples) the lags for xcorr

% n_fft: number of frequency bands
[n_bands, n_smp] = size(X);

% Calculate the bias of the cross-correlation ('unbias') 
bias = xcorr(ones(1,n_smp), ones(1,n_smp), maxlag);



%% AUTO-CORRELATION vs. CROSS-CORRELATION
if isvector(Y)
    % Stimulus\Stimulus AUTO-CORRELATION Case 
    
    % # of cross-correlations == # of frequency bands
    n_xcorr = n_bands;
    
    y_idx = ones(n_xcorr, 1);
    x_idx = 1:n_bands;
    
else 
    % Stimulus\Response CROSS-CORRELATION 

    % Indices to the correlations
    [y_mesh, x_mesh] = meshgrid(1:n_bands);

    % The correlation matrix is conjugate-symmetric, so we need to calculate
    % only its lower side + it diagonal
    y_idx = nonzeros( tril(y_mesh) );
    x_idx = nonzeros( tril(x_mesh) );

    % # of correlations (# of upper entries + diagonal)
    n_xcorr = n_bands*(n_bands+1)/2;
    
end



%% CORRELATION MATRIX
switch xcorr_type
    case 1
        % Option #1 (FFT); slow

        % First calculate all Fourier transforms, and then perform the
        % assignment
        Xw = fft(X, 2*n_smp+1, 2);     % 2*n_smp: avoid cyclic aliasing
        Yw = fft(Y, 2*n_smp+1, 2);     % 2*n_smp: avoid cyclic aliasing
        
        % See Theunissen et al., 2000, <X(w),Y*(w)>
        Cab = ifft( Xw(x_idx,:) .* conj(Yw(y_idx,:)), [], 2, 'symmetric');  

        % Apply the window
        Cx = [Cab(:, end-n_win+1:end), Cab(:, 1:n_win+1)];  % fftshift on maxlags

        % UNBIAS the correlations; 
        if remove_bias
            % REMOVE the bias (it's like doing xcorr and setting scaleopt to
            % 'unbiased')
            Cx = Cx./bias;
        end
    
    
    case 2
        % Option #2 (regular XCORR); slowest
        if remove_bias
            scaleopt = 'unbiased'; 

            % In this case there is no bias for the XCORR; normalization is
            % performed by division, so BIAS == 1 <==> no bias.
            bias = ones(size(bias));  
        else
            scaleopt = 'none'; 
        end
        
        Cx = zeros(n_xcorr, 2*maxlag+1);
        for k = 1:n_xcorr
            I1 = x_idx(k);
            I2 = y_idx(k);
            Cx(k,:) = xcorr(X(I1,:), Y(I2,:), scaleopt, maxlag );
        end


    case 3
        % Option #3: (regular XCORR) Get the correlation entries along the "long" 
        % dimension, so the for loop is "shorter" (i.e., along MAXLAG). 
        % My reference for this code is taken from the STRFLAB project
        % (see df_small_autocorr4.m).                
        
        Cx = zeros(n_xcorr, 2*maxlag+1);  
        
        % *** USE GPU ***
        if strcmpi('gpuArray',variable_type)
            Cx = gpuArray(Cx);
            X  = gpuArray(X);
            Y = gpuArray(Y);
        end
        
        idx_get_tril = 1 == tril(ones(size(X,1), size(Y,1)));        
        for tid = -n_win:n_win
            v = (max(1,tid+1)):(min(n_smp,n_smp+tid));
            w = v - tid;
            temp = X(:,v) * Y(:,w)';
            Cx(:,n_win + tid + 1) = temp(idx_get_tril); 
        end
        
        % *** USE GPU ***
        if strcmpi('gpuArray',variable_type)
            Cx = gather(Cx);
        end
        
        % UNBIAS the correlations; 
        if remove_bias
            % REMOVE the bias (it's like doing xcorr and setting scaleopt to
            % 'unbiased')
            Cx = Cx./bias;
        end
        
        
    otherwise
        error('--> [autocorr_mtx.m]: unrecognized XCORR type (xcorr_type: %d)!!\n', xcorr_type);
        
end
    









