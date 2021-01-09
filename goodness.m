function gof = goodness(X, Y, type_out)
%
%   [CC, mse, nmse] = goodness(x, y, type_out)
%

if 3 > nargin
    type_out = '';
end

% if 4 > nargin
%     apply_whitening = true;
% end


%% Set the compared inputs (matrices or vectors) into columns vectors
x = X(:);
y = Y(:);

% Remove the mean
x = x - mean(x);
y = y - mean(y);

%% CC
R = corrcoef( x, y );
gof.CC = R(1,2);
if strcmpi('CC', type_out)
    gof = gof.CC;
    return;
end


%% MSE & NMSE
% %{
SS   = @(a) mean(abs( a - mean(a) ).^2);
% mse1 = @(a) mean(abs( a ).^2);
mse2 = @(a,b) mean(abs( a - b ).^2);

gof.mse = mse2(x, y);
gof.nmse = gof.mse/SS(x);
% gof.nmse= gof.mse/sqrt(mse1(x).*mse1(y));
% gof.nmse= gof.mse/sqrt( max(x.^2) - min(x.^2) );
%}


%{
rmsd = @(a,b) sqrt(mean(( a - b ).^2));   % root-mean-square deviation
nrmsd= @(a,b) rmsd(a,b) ./ ( max(abs(a)) - min(abs(a)) );

gof.mse = rmsd(x, y);
gof.nmse= nrmsd(x, y);
if strcmpi('mse', type_out)
    gof = gof.mse;
    return;
elseif strcmpi('nmse', type_out)
    gof = gof.nmse;
    return;    
end
%}



%% Mean correlation coefficients over frequencies
% Noise-invariant Neurons in the Avian Auditory Cortex:
% Hearing the Song in Noise; Moore, Theunissen et al., 2013.
X_minus_mu = X - mean(X,2);
Y_minus_mu = Y - mean(Y,2);
% X_minus_mu = zca(X')';
% Y_minus_mu = zca(Y')';
num = mean(X_minus_mu .* Y_minus_mu, 2);
den = sqrt( var(X_minus_mu,[],2) .* var(Y_minus_mu,[],2) );     % den ~ 1.0 due to the zca()
gof.CCf = num ./ den;
if strcmpi('CCf', type_out)
    gof = gof.CCf;
    return;
end




