function [D, dict_st, Dimg] = compute_bases(Sft, n_lags, varargin)
%
%   function [D, dict_st, Dimg] = compute_bases(Sft, n_lags, ...)
%
%

p = inputParser;

% === Required inputs ===
addRequired(p, 'Sft', @isnumeric);              % (n_bands x n_time) spectrogram power
% addRequired(p, 'n_bands', @isnumeric);            
addRequired(p, 'n_lags', @isnumeric);  	

% === Optional inputs ===
addParameter(p, 'dict_type', 'FASTICA', @ischar);    % {'FASTICA', 'RICA', 'PCA', 'DL', 'NNMF'}
addParameter(p, 'n_atoms', 9, @isnumeric);           
addParameter(p, 'f_bpf', [], @isnumeric);        % bandpass filter
addParameter(p, 'whitening', true, @(x) islogical(x) || isnumeric(x));  

% DL:
addParameter(p, 'lambda', 0.01, @isnumeric);     % (1x1) sparse coding coefficient    

% For the plots\debug
addParameter(p, 'verbose', 0, @isnumeric);           
addParameter(p, 'fignum', [], @isnumeric);           

% === Parse ===
parse(p, Sft, n_lags, varargin{:});
pars = p.Results;


% Whitening (sphering)
% zca = @(x) (x-mean(x))./std(x);

% Initiate the matrix of the dictionary
D = [];

% number of frequency channels
n_bands = size(Sft, 1);

if isempty(pars.fignum)
    pars.fignum = 0;
end


%% INit. the band-pass filter (if used)
% Apply an octave-band filter, if pars.f_bpf frequency is supplied;
bpf = [];
if ~isempty(pars.f_bpf)
    % Octave filter
    bpf.bw = 1;                 % octave
    bpf.filter_order  = 4;      % (1x) filter's order

    fs = 4 * max(pars.f_bpf);       % (Hz) times 4 for the bandpass-octave filter
    fc = pars.f_bpf;  % (Hz)
    [bpf.B, bpf.A] = octave.design_filter(fc, fs, bw, filter_order);
    [Hk, f_] = freqz(bpf.B, bpf.A, pars.f_bpf, fs);
    assert(all( f_ == pars.f_bpf ));
    assert( ~isnan(sum(abs(Hk))) );
    Sft = Sft .* abs(Hk);
    aux.vprint(pars.verbose, '[compute_bases.m]: Using an octave-bandpass filtering (fc: %g Hz)\n', fc);                        
end



%% Cut chunks out of the spectrogram 
X = im2col(Sft, [n_bands, n_lags], 'sliding');

% % Remove from each patch (column) its mean
% X = X - repmat(mean(X), [size(X,1) 1]);
% 
% % Normalize each patch by its variance
% X = X ./ repmat(eps + sqrt(sum(X.^2)),[size(X,1) 1]);     

% # of features\bases in the dictionary
pars.n_atoms = min(size(X,2), pars.n_atoms);      
pars.dict_type = upper(pars.dict_type);
dict_st = pars;



%%
% Remove from each patch (column) its mean
X = X - repmat(mean(X), [size(X,1) 1]);

% Normalize each patch by its variance
X = X ./ repmat(eps + sqrt(sum(X.^2)),[size(X,1) 1]);    

switch pars.dict_type
    case 'RICA'
        ica_model = rica(X', pars.n_atoms, 'IterationLimit', 1000, 'Lambda', dict_st.lambda);
        D = ica_model.TransformWeights;   % (img x features)
        
        dict_st.ica_model = ica_model;
        
    case 'FASTICA'   
        % Remove from each patch (column) its mean
        %X = X - repmat(mean(X), [size(X,1) 1]);
        %
        % Normalize each patch by its variance
        %X = X ./ repmat(eps + sqrt(sum(X.^2)),[size(X,1) 1]);     

        dict_st.approach    = 'symm';    % {'defl', 'symm'}
        dict_st.only        = 'all';
        dict_st.maxFinetune = 500;
        dict_st.n_atoms  = pars.n_atoms;
        verbose     = 'Off';

        [~, D, ~] = fastica.fastica(X, ... (:, 1:max_bases),...
            ... 'lastEig', 10, ...
            'approach', dict_st.approach, ... {'defl', 'symm'}
            'only', dict_st.only, ...
            'maxFinetune', dict_st.maxFinetune, ...
            'numOfIC', dict_st.n_atoms, ...
            'verbose', verbose );



    case 'PCA'
        % Remove the mean from each patch (column) 
        %X = X - repmat(mean(X), [size(X,1) 1]);

        % Normalize each patch (column) by the variance
        %X = X ./ repmat(eps + sqrt(sum(X.^2)),[size(X,1) 1]);     

        [D,S,~] = svds(X, pars.n_atoms);        
        dict_st.S = diag(S);    % save the singular values


    case 'DL'   % dictionary learning
        % Remove from each patch (column) its mean
        %X = X - repmat(mean(X), [size(X,1) 1]);

        % Normalize each patch by its variance
        %X = X ./ repmat(eps + sqrt(sum(X.^2)),[size(X,1) 1]);    

        dict_st.K           = pars.n_atoms;	% learns a dictionary with K elements
        %dict_st.lambda      = pars.lambda;
        dict_st.numThreads  = -1;           % number of threads
        %dict_st.batchsize   = 400;          % default is 512
        dict_st.iter        = 1000;         % let us see what happens after 1000 iterations.
        %pars.verbose        = false;

        aux.vprint(pars.verbose, '--> Starting dictionary learning (DL)...\n');
        D = mexTrainDL(X, dict_st);
        aux.vprint(pars.verbose, '--> Finished\n');



    case 'NNMF'   % non-negative matrix factorization
        % k must be a positive integer no larger than the number of rows or columns in A
        new_n_features = min(pars.n_atoms, size(X,1));
        new_n_features = min(new_n_features, size(X,2));
        if pars.n_atoms ~= new_n_features
            warning('--> [compute_bases.m]: n_atoms (%d) was changed to (%d)', pars.n_atoms, new_n_features);
            pars.n_atoms = new_n_features;
        end
        
        try
            [D, ~] = nnmf(X, pars.n_atoms);
        catch
            D = [];
        end

    otherwise
        error('-> Unrecognized dict_type!!');
end


empty_columns = 0==sum(D);
if  0 < nnz(empty_columns)
    warning('---> [compute_bases.m]: There are empty bases (atoms) in this dictionary!!');
    fprintf('---> Empty dictionary columns are removed from the dictionary!\n');
    D = D(:, ~empty_columns);
    
    % update number of features/atoms
    dict_st.n_atoms = size(D, 2);
end




% 'DEBUG'
%{
'*** DEBUG ***'
figure(1);
clf
subplot(1,2,1);
imagesc(reshape(D(:,1), n_bands, n_lags));
%figure(2);
%clf
subplot(1,2,2);
loglog(pars.f_bpf, abs(Hk)); grid on    
title(sprintf('k: %d', k));
drawnow;
pause(0.01);
%}

% Whitening each column (atom) in the dictionary
if pars.whitening
    D = zca(D);
end




%% Reshape the columns back into 2D filters; 
% - Dimg is a tensor of size (num of atoms, num frequency bands, num of temporal lags)
dict_st.atom_dim1 = n_bands;
dict_st.atom_dim2 = n_lags;

if 2 < nargout
    % # of filters (bases)
    pars.n_atoms = size(D, 2);

    Dimg = nan(pars.n_atoms, n_bands, n_lags);
    for k = 1:pars.n_atoms
        Dimg(k,:,:) = reshape(D(:,k), n_bands, n_lags);
    end
    
end




%% Print & plot stuff
if pars.verbose
    aux.printstruct(dict_st);
end

if ~pars.fignum
    return; 
end
figure(pars.fignum);
clf;

% plot_dic_patches( A(:,1:(9*10)), n_bands, n_lags );
[ax, subs_x, ~] = plot_dic_patches( D, n_bands, n_lags );
title(ax(ceil(subs_x/2)), sprintf('Selected Filters (%s)', pars.dict_type));
















