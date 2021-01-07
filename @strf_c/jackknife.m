function [strf, JN] = jackknife(obj, X_train, r_train, varargin)
%
%   function JN = jackknife(X_train, r_train, varargin)
%


%% Set the input
p = inputParser;


% === Required inputs ===
addRequired(p, 'X_train', @isnumeric);      % (n_bands x n_time) spectrogram power
addRequired(p, 'r_train', @isvector);       % (n time x 1) psth of the response


% === Optional inputs ===
% Algorithm type to calculate the STRF 
%addOptional(p, 'algo_type', obj.algo_type, @isstr);      % % {'ASD', 'regression'}             

% addOptional(p, 'lags_ms', obj.lags_ms, @isscalar);      % # of chunks for the jackknife   
addOptional(p, 'n_splits', 11, @isscalar);      % # of chunks for the jackknife   
addOptional(p, 'smooth_response_flag', obj.smooth_response_flag, @isscalar);      % # of chunks for the jackknife   

% Tolerance for the ridge regression. If tol is an array, the function will 
%   return an array (s) of STRFs that relates to each of the M values.
% addOptional(p, 'gamma', obj.ridge.gamma, @isnumeric);      % (Mx1)       
addOptional(p, 'gamma_4_JK', obj.ridge.gamma_4_JK, @isnumeric);      % (Mx1)       
 
% Shrinkage parameter for the shrinkage filter
% - Stephen David et al., 2004, Natural Stimulus Statistics Alter the 
%   Receptive Field structure of V1 neurons
addOptional(p, 'theta_4_JK', obj.theta_4_JK, @isnumeric);      % (Mx1)       

addOptional(p, 'xcorr_type', obj.xcorr_type, @isnumeric);    % (1x1) correlation method       

% If 1 (true) then output a causal STRF(s), that is, all STRFs are zero for
%   t<0 ms. Otherwise, if 0, return the "non-causal" STRF(s).
addOptional(p, 'iscausal', obj.iscausal, @isnumeric);      % (1x1)       


addOptional(p, 'verbose', 0, @isnumeric);       % (1x1) write things to the command line?
addOptional(p, 'fignum', [], @isnumeric);       % (1x1) if not empty, plot stuff (for debug)

parse(p, X_train, r_train, varargin{:});
pars = p.Results;


% Save changes into the object, if needed
%obj.algo_type           = pars.algo_type;
obj.iscausal            = pars.iscausal;  
n_splits                = pars.n_splits;     
obj.ridge.gamma_4_JK    = pars.gamma_4_JK;     
obj.smooth_response_flag= pars.smooth_response_flag;
gamma_4_JK              = obj.ridge.gamma_4_JK;           % for a short access

% Stephen David et al., 2004, Natural Stimulus Statistics Alter the 
%   Receptive Field structure of V1 neurons
obj.theta_4_JK = pars.theta_4_JK;       
obj.xcorr_type = pars.xcorr_type;   


% STRF parameters
n_lags     = obj.n_lags;
n_bands    = size(X_train, 1);
n_gamma    = numel(gamma_4_JK);          
n_theta    = numel(obj.theta_4_JK);

fignum     = pars.fignum;       
verbose    = pars.verbose;

obj.X_train = X_train;
obj.r_train = r_train;
assert(strcmpi('regression', obj.algo_type),...
    '--> [strf_c/jackknife.m]: currently you can''t do jackknife with ASD; You need to change ALGO_TYPE back to ''regression''!');

if 1 == verbose 
    fprintf('\n-> Starting jacjknife...\n');
    fprintf('  -----------------------\n');
    fprintf('--> iscausal    : %d\n', obj.iscausal);    
    fprintf('--> # ridge tol : %d\n', n_gamma);
    fprintf('--> ridge to try: [%s]\n', num2str( gamma_4_JK, ' %g'));    
    fprintf('--> # shrinkage : %d\n', n_theta);
    fprintf('--> n_lags      : %d\n', obj.n_lags);
    fprintf('--> xcorr_type  : %d\n', obj.xcorr_type);
end



%% Cost Function
if pars.smooth_response_flag
    aux.cprintf('r', '--> [jackknife.m]: smooth_response_flag is ON!\n');
    
    % Add a window to smooth the response(especially for the PSTH)
    cost.win_size       = floor(obj.n_lags/2);                           % (samples)
    cost.win_fun        = @(n) gausswin(n);             % gausswin(n)\ hann(n)
    cost.win            = cost.win_fun( cost.win_size );
    cost.win            = cost.win/sum(cost.win);       % normalize
else
    cost.win            = [];     % no smoothing window
end

cost.distance       = 'correlation';  % {'euclidean', 'seuclidean', 'cosine', 'correlation'}
cost.distParameter  = [];
cost.fun = @(x,y) pdist([x(:), y(:)].', cost.distance, cost.distParameter);



%% Split the data
% Split the data into (approximately) equal exclusive sets for the jackkinfe
dim_x           = 2;
nt              = size(X_train, dim_x);  % samples along the time domain
splits.idx      = fix((0:(nt-1))/ceil(nt/n_splits)) + 1;
splits.n_splits = n_splits;
splits.grps     = 1:n_splits;



%% Jackknife structure
% JN.tol      = gamma_4_JK;
JN.gamma_4_JK= obj.ridge.gamma_4_JK;
JN.theta_4_JK= obj.theta_4_JK;
JN.splits   = splits;
JN.n_splits = n_splits;
JN.strf_info= 'Dims: (n_bands, n_lags, n_splits, n_tol)';
JN.strf     = nan(n_bands, n_lags, n_splits, n_gamma);  % mean(STRFs) across the JN set, given the ridge's tolerance
JN.mean     = nan(n_bands, n_lags, n_gamma);  % mean(STRFs) across the JN set, given the ridge's tolerance
JN.SE       = nan(n_bands, n_lags, n_gamma); % std(STRFs) across the JN set, given the ridge's tolerance
JN.avg_psth = zeros(1, n_gamma);  % 
JN.sparse_strf = nan(n_bands, n_lags, n_theta, n_gamma);  % mean(STRFs) across the JN set, given the ridge's tolerance

JN.error_info = 'Dims: (n_gamma x n_tol)';
JN.error = zeros(n_theta, n_gamma);  

JN.cost = cost;



%%
aux.vprint(verbose, '\n');
for mm = 1:n_gamma
    aux.vprint(verbose, '---> [jackknife.m]: tolerance: %f (%d out of %d)\n',...
        gamma_4_JK(mm), mm, n_gamma);

    for nn = 1:n_splits     % # of splits\chunks
        % Select all except the nn chunk
        test_number = nn;
        train_idx = splits.idx ~= test_number;
        Xn = X_train(:, train_idx);
        yn = r_train(train_idx);
        
        % Calculate STRF of 
        strf_mn = obj.calc_strf(Xn, yn, 'gamma', gamma_4_JK(mm), 'algo_type', obj.algo_type);
        
        % Save the STRF without the nn'th part
        JN.strf(:,:,nn,mm) = strf_mn;
        
        % Add back the average of the PSTH (without the bb'th part)
        %JN.avg_psth(mm) = JN.avg_psth(mm) + strf_mn_st.avg.psth;
        JN.avg_psth(mm) = JN.avg_psth(mm) +  obj.bias.psth;        
    end
    
    % Average over the sum of the PSTHs to get the mean value
    JN.avg_psth(mm) = JN.avg_psth(mm)/JN.n_splits;
    
    % mean(STRFs): the mean STRF across the JN set given the ridge's tolerance
    JN.mean(:,:,mm) = mean( JN.strf(:,:,:,mm) , 3 );
    
    % std(STRFs): the mean STRF across the JN set given the ridge's tolerance
    JN.SE(:,:,mm) = std( JN.strf(:,:,:,mm), [], 3 );
    
    % Run over the optimal shrinkage parameter
    for gg = 1:n_theta
        JN.sparse_strf(:,:,gg,mm) = shrinkage_fun(...
            JN.mean(:,:,mm),...     mean STRF for a given tolerance (tolerance: ridge regression ocef.)
            JN.SE(:,:,mm),...       SE (std) of the mean STRF for a given tolerance
            obj.theta_4_JK(gg)...   the 'shrinkage', or 'sparse' (aka STRFLAB) coefficient
            );

        for kk = 1:JN.n_splits
            test_idx = splits.idx == kk;
            Xn = X_train(:, test_idx);
            yn = r_train(test_idx);
            
            % Optimal linear response estimation
            rest_gn = obj.predict(Xn, JN.sparse_strf(:,:,gg,mm),...
                'avg', JN.avg_psth(mm), ...
                'normalize', max(yn) );

            % !!! EXPERIMENTAL !!!
            if ~isempty(cost.win)
                yn_filt = filtfilt(cost.win, 1, yn);
            else
                yn_filt = yn;   % no filtering
            end
            JN.error(gg, mm) = JN.error(gg, mm) + cost.fun(yn_filt, rest_gn);
        end
        
        % Average the cost over the number of jackknife cases
        JN.error(gg, mm) = JN.error(gg, mm)/JN.n_splits;
        
    end   
    
end

% Save all results\outputs
[JN.best.error, best_idx]       = min(JN.error(:));
[JN.best.sub(1), JN.best.sub(2)]= ind2sub(size(JN.error), best_idx);
JN.best.gamma                   = gamma_4_JK(JN.best.sub(2));
JN.best.theta                   = obj.theta_4_JK(JN.best.sub(1));
JN.best.avg_psth                = JN.avg_psth( JN.best.sub(2) );

% The best STRF
% JN.best.strf = JN.sparse_strf(:,:,JN.best.sub(1),JN.best.sub(2));
strf    = JN.sparse_strf(:,:,JN.best.sub(1),JN.best.sub(2));
obj.strf= strf;



%% Plot
if isempty(fignum), return; end

figure(fignum);
clf;
% subplot(1,2,1);
imagesc(log10(gamma_4_JK), log10(obj.theta_4_JK), JN.error );
hold on
plot(log10(JN.best.gamma), log10(JN.best.theta), 'xr');
hold off
colorbar;
set(gca, 'YDir', 'normal');
ylabel('$\theta$ Shrinkage Coef. $(\times10^{-3})$');
xlabel('Ridge Tolerance $(\times10^{-3})$');
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%.2f',x), 1e3*obj.theta_4_JK, 'UniformOutput', 0));
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.2f',x), 1e3*sort(obj.ridge.gamma_4_JK), 'UniformOutput', 0));
title(sprintf('Discrepancy (best shrinkage coef: %g, best ridge tol: %g)',...
    obj.theta_4_JK(JN.best.sub(1)), gamma_4_JK(JN.best.sub(2))));

figure(1+fignum);
clf;
ax = obj.plot_strf;
title(ax(1), 'Best STRF');











