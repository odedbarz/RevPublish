function G = jackknife(obj, X_train, r_train, varargin)
%            
% function G = jackknife(obj, X_train, r_train, [n_splits], [G], [lags], [bias])
% 


%% Set the input data
p = inputParser;

% Required
% addRequired(p, 'obj', @isobject);            
addRequired(p, 'X_train', @isnumeric);            
addRequired(p, 'r_train', @isnumeric);            

% Optional
addOptional(p, 'lags_ms', obj.lags_ms, @isnumeric);         % (ms)
addOptional(p, 'f', obj.f, @isnumeric);                     % (Hz) frequency axis     
addOptional(p, 'n_splits', obj.n_splits, @isnumeric); 
addOptional(p, 'algo_type', obj.algo_type, @isstr); 
addOptional(p, 'iscausal', obj.iscausal, @(x) isnumeric(x) || islogical(x)); 
addOptional(p, 'gamma_4_JK', obj.ridge.gamma_4_JK, @isnumeric);           
addOptional(p, 'verbose', 0, @isnumeric);           
addOptional(p, 'fignum', [], @isnumeric);           


parse(p, X_train, r_train, varargin{:});
pars = p.Results;


%%
% Update the object, if needed
obj.lags_ms         = pars.lags_ms;     
obj.iscausal        = pars.iscausal;
obj.n_splits        = pars.n_splits;     
obj.algo_type       = pars.algo_type;     

% In case of algo_type == 'svd' there is no need to run various gamma
% values!
if strcmpi('svd', obj.algo_type) && length(obj.ridge.gamma_4_JK) > 1
    aux.vprint(pars.verbose, '--> [reconstruct_c/jackknife.m]: when using ALGO_TYPE == ''svd'' use only one gamma value!\n');
    obj.ridge.gamma_4_JK = obj.ridge.gamma;
else
    obj.ridge.gamma_4_JK = pars.gamma_4_JK;    
end


obj.X_train         = X_train;
obj.r_train         = r_train;

verbose    = pars.verbose;                
fignum     = pars.fignum;       

% For the current function
n_splits   = pars.n_splits;
gamma_4_JK = obj.ridge.gamma_4_JK;   % a list for the ridge regression jackknife
assert(1 < n_splits );

% Split the data into (approximately) equal exclusive sets for the jackkinfe
nt              = size(X_train, 2);  % samples along the time domain
splits.idx      = fix((0:(nt-1))/ceil(nt/n_splits)) + 1;
splits.n_splits = n_splits;
splits.grps     = 1:n_splits;


% Dimensions for the G matrix
n_bands= size(X_train, 1);     % # of frequency bands
n_lags = obj.n_lags;
dim1_G = n_lags * obj.n_neurons;



%% BIAS -- TRAINING -- REMOVE
% if 1 == obj.bias_choice
%     obj.bias.remove = 1;
% end

% Keep the bias for the fitting 
obj.bias.train.Sft  = mean(X_train,2); 	% spectrogram's bias 
obj.bias.train.resp = mean(r_train);      % PSTH's bias

n_gammas= length(gamma_4_JK);
G_mean  = nan(dim1_G, n_bands, n_gammas);
G_sd    = nan(dim1_G, n_bands, n_gammas);
scores  = nan(1, n_gammas);
Crr     = nan(dim1_G, dim1_G, n_gammas);
Crr_inv = nan(dim1_G, dim1_G, n_gammas);
Crs     = nan(dim1_G, n_bands, n_gammas);



%% Start jackknife
aux.vprint(verbose, '--> [JK.m]: starting the jackknife for loop over the GAMMAs...\n');
for m = 1:n_gammas
    G_k     = nan(dim1_G, n_bands, n_splits);
    Crr_k   = nan(dim1_G, dim1_G, n_splits);
    Crr_inv_k= nan(dim1_G, dim1_G, n_splits);
    Crs_k   = nan(dim1_G, n_bands, n_splits);
    score_k = nan(1, n_splits);
        
    % Construct the reconstruction object
    obj_m = reconstruct_c(obj.binwidth,...
        'f', obj.f, ...
        'lags_ms', obj.lags_ms,...
        'algo_type', obj.algo_type, ...
        'iscausal', obj.iscausal, ...
        'inv_type', obj.ridge.inv_type, ...
        'gamma', gamma_4_JK(m),...
        'remove_bias', obj.bias.remove,...
        'gpu_flag', obj.gpu_flag, ...
        'xcorr_type', obj.ridge.xcorr_type); 
    aux.vprint(verbose, '--> [JK.m]: gamma: %g\n', obj.ridge.gamma);

    for k = 1:n_splits
        aux.vprint(verbose, '--> [JK.m]: split # %d (%d)\n', k, n_splits);

        test_grp_k= splits.grps(k);
        X_train_k = X_train(:, splits.idx ~= test_grp_k);
        r_train_k = r_train(splits.idx ~= test_grp_k, :);

        % Fit the model (find the "inverse-STFT" filters G)
        obj_m.calc_filter(X_train_k, r_train_k);
        
        assert(~isempty(obj_m.G), '--> [reconstruct_c/jackknife.m]: The reconstruction filters is invalid!');
        G_k(:,:,k) = obj_m.G;
        
        % Only the ALGO_TYPE == 'regression' calculates these matrices
        if ~isempty(obj_m.Crr)
            Crr_k(:,:,k)    = obj_m.Crr;
            Crr_inv_k(:,:,k)= obj_m.Crr_inv;
            Crs_k(:,:,k)    = obj_m.Crs;
        end
        
        % Predict the spectrogram
        rk_test_k = r_train(splits.idx == test_grp_k, :);
        X_test_k  = X_train(:, splits.idx == test_grp_k);
        Xk_est    = obj_m.predict(rk_test_k);
        
        score_k(k) = goodness(X_test_k, Xk_est, 'CC');
                
    end
    
    % Average over all the splits for a given gamma
    G_mean(:,:,m) = mean(G_k, 3);
    G_sd(:,:,m)   = std(G_k, [], 3);
    scores(m)     = median(score_k);
    Crr(:,:,m)    = median(Crr_k, 3);
    Crr_inv(:,:,m)= median(Crr_inv_k, 3);
    Crs(:,:,m)    = median(Crs_k, 3);
    
end
aux.vprint(verbose, '--> [JK.m]: Finished the jackknife for loops\n');



%% Save the best setup into the object
[best_score, best_score_idx] = max(scores);     % get the G with the best score (using jackkife)
obj.G                   = G_mean(:,:,best_score_idx);
obj.G_sd                = G_sd(:,:,best_score_idx);
obj.Crr                 = Crr(:,:,best_score_idx);
obj.Crr_inv             = Crr_inv(:,:,best_score_idx);
obj.Crs                 = Crs(:,:,best_score_idx);
obj.ridge.scores        = scores;
obj.ridge.best_gamma    = gamma_4_JK(best_score_idx);
obj.ridge.best_score    = best_score;
obj.ridge.best_score_idx= best_score_idx;



%% Set the output, if required
if 0 < nargout
    G = obj.G;     
end



%% Plot the reconstructed spectrogram
if isempty(fignum) || fignum<=0, return; end

% (Hz or samples)
if isempty(obj.f)
    f = 1:obj.n_bands;  % (samples)
    ylabel_str = 'Frequency (idx)';
else
    f = 1e-3*obj.f;          % (Hz)
    ylabel_str = 'Frequency (kHz)';    
end

figure(fignum);
clf;

% Plot 1 example:
ax = subplot(1,2,1);
surf_h = surf(obj.binwidth*obj.lags, f, obj.G(1:n_lags,:)');
view(2);
axis tight
set(surf_h, 'EdgeColor', 'none');
colorbar;
title('$mean(G_1)$');
xlabel('Lags (ms)');
ylabel( ylabel_str );
colormap jet

ax(2) = subplot(1,2,2);
surf_h = surf(obj.binwidth*obj.lags, f, obj.G_sd(1:n_lags,:)');
view(2);
axis tight
set(surf_h, 'EdgeColor', 'none');
colorbar;
title('$SD(G_1)$');
xlabel('Lags (ms)');
% caxis(ax(2), caxis(ax(1)));
colormap jet

figure(2+fignum);
clf;
semilogx(obj.ridge.gamma_4_JK, scores, '.:')
title(sprintf('Jackknife Scores $(\\gamma_{opt}: %g)$', obj.ridge.best_gamma));
xlabel('$\gamma$');
ylabel('score');
hold on
plot(obj.ridge.gamma_4_JK(best_score_idx), obj.ridge.best_score, 'rs');
hold off





