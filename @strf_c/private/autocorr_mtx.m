function [Cx, bias] = autocorr_mtx(Sx, n_win, xcorr_type, remove_bias)
%
% function [Cx, bias] = corr_mtx(Sx, n_win, xcorr_type)
%

if 3 > nargin || isempty(xcorr_type)
    xcorr_type = 1;
end

if 4 > nargin || isempty(remove_bias)
    remove_bias = false;     % remove the bias (unbias) from the cross-correlation
end

%% Set the input variables
% Correlation parameters
maxlag = n_win;      	% (samples) the lags for xcorr

% n_fft: number of frequency bands
[n_fft, n_smp] = size(Sx);

% # of correlations (# of upper entries + diagonal)
n_xcorr = n_fft*(n_fft+1)/2;

% Indices to the correlations
[x_mesh, y_mesh] = meshgrid(1:n_fft);

% The correlation matrix is conjugate-symmetric, so we need to calculate
% only its lower side + it diagonal
x_idx = nonzeros( tril(x_mesh) );
y_idx = nonzeros( tril(y_mesh) );

% % Remove stimulus means along the y-axis (frequency bands)
% Sft = Sft - mean(Sft,2);

% Calculate the bias of the cross-correlation ('unbias') 
bias = xcorr(ones(1,n_smp), ones(1,n_smp), maxlag);

%% CORRELATION MATRIX
% %{
switch xcorr_type
    case 1
        % Option #1 (FFT); slow

        % First calculate all Fourier transforms, and then perform the
        % assignment
        Sw = fft(Sx, 2*n_smp+1, 2);     % 2*n_smp: avoid cyclic aliasing
        
        % See Theunissen et al., 2000, <S*(w),S(w)>
        Cab = ifft( Sw(y_idx,:) .* conj(Sw(x_idx,:)), [], 2, 'symmetric');  

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
            I1 = y_idx(k);
            I2 = x_idx(k);
            Cx(k,:) = xcorr(Sx(I1,:), Sx(I2,:), scaleopt, maxlag );
        end


    case 3
        % Option #3: (regular XCORR) Get the correlation entries along the "long" 
        % dimension, so the for loop is "shorter" (i.e., along MAXLAG). 
        % My reference for this code is taken from the STRFLAB priject
        % (see df_small_autocorr4.m).        
        Cx = zeros(n_xcorr, 2*maxlag+1);        
        Sx_ = Sx';
        for tid = -n_win:n_win
            v = (max(1,tid+1)):(min(n_smp,n_smp+tid));
            w = v - tid;
            temp = Sx(:,v)*Sx_(w,:);
            Cx(:,n_win + tid + 1) = nonzeros(tril(temp,0));
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
    









