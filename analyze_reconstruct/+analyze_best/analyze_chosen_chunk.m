% analyze_chosen_chunk.m
%
% Analyze a chosen chunk.


% Train: orthogonalize the measurements and get the most effective set of units
train_drr   = 3     % drr.dry;
test_drr    = 4     % drr.dry;
test_speaker= 5;    % chunk\speaker 
n_units     = 20;   % Choose # of units to use for the reconstruction

aux.cprintf('string', '--> train_drr   : %d\t (%s)\n', train_drr, drr.labels{train_drr});
aux.cprintf('string', '--> test_drr    : %d\t (%s)\n', test_drr, drr.labels{test_drr});
aux.cprintf('string', '--> test_speaker: %d\t (%s)\n', test_speaker, tbl_metadata.fn{test_speaker});

% sort_type = 'ascend';
% sort_type = 'randperm';
sort_type = 'orthogonal';
% sort_type = 'invariant-strf';

H_train = squeeze( H_unit(:, train_drr, :) );


% Get an optimal set
switch sort_type
    case 'orthogonal'         
        '!!! Over-write !!!'
        aux.cprintf('string', '--> sort_type: %s, !!! Over-write train_drr to select units !!!\n', sort_type);
        H_ = squeeze( H_unit(:, drr.dry, :) );   % Over-write
        [slc.optimal_sorted, P, sv] = find_best_unit_set(H_, 'n_svd', 36);
        
        % transpose to get the neurons for each orthogonal vector
        slc.optimal_sorted = unique(slc.optimal_sorted', 'stable');

    case 'randperm'
        slc.optimal_sorted = randperm(size(H_unit,3))';     % a random permutation
        aux.cprintf('r', '--> sort_type: RANDOM\n');
        
    case 'ascend'
        slc.optimal_sorted = (1:size(H_unit,3))';     % a ascending order
        aux.cprintf('r', '--> sort_type: ASCEND\n');
        
    case 'invariant-strf'
        dummy = load('../.data/STRFs/invariant_bf_list.mat');
        slc.optimal_sorted = find(dummy.inv_bf);
        
    otherwise
        error('--> Unrecognized sort type!');
end

% Use these units
%slc.optimal_sorted = slc.optimal_sorted(:)';
slc.optimal_units_used = slc.optimal_sorted(1:n_units, :);
aux.cprintf('string', '--> # of units: %d\n', n_units);

% % Training responses
% %y1 = squeeze( H_unit(:, train_drr, slc.optimal_units_used) );
% y1 = H_train(:, slc.optimal_units_used);
% 
% % Training spectrogram
% X1 = spec_st.Sft{train_drr};

% Plot the projections 
% %{
if 1 == verbose
    fprintf('--> slc.optimal_units_used: ');
    disp(slc.optimal_units_used(:)');
end

if ~isempty(fignum) && strcmpi(sort_type, 'orthogonal')
    figure(fignum);
    clf;
    plot(P(:,1), '.-');
    xlabel('Unit Number');
    ylabel('$Y_1^T\cdot Z$');
    title('Projections over Orthogonal Axes');
end
%}



%% Train the model using the "best" chosen units
X1 = spec_st.Sft{train_drr};                % trained spectrogram
y1 = H_train(:, slc.optimal_units_used);    % trained responses

[X_train, X_test0, y_train, y_test0, ~] = train_test_split(X1, y1, ...
    'split_time_idx', split_time_idx,...
    'test_grp', test_speaker ...
);

% Initialize the reconstruction object
obj = reconstruct_c(binwidth,...
    'f', spec_st.f,...
    'iscausal', iscausal, ...
    'algo_type', algo_type, ...
    'lags_ms', lags_ms ); 

% Fit the model
obj.fit(X_train, y_train,...
    'jk_flag', jk_flag,...
    'n_splits', n_splits, ...
    'fignum', []);



%% Test 
X2 = spec_st.Sft{test_drr};                                     % Testing spectrogram
y2 = squeeze( H_unit(:, test_drr, slc.optimal_units_used) );    % Testing responses

[~, X_test, ~, y_test, ~] = train_test_split(X2, y2, ...
    'split_time_idx', split_time_idx,...
    'test_grp', test_speaker );    


% PREDICTION; Predict the spectrogram using the selected responses
X_est = obj.predict(y_test);    
gof = goodness(X_test, X_est);



%% Plot stuff
if verbose
    fprintf('\n');
    if ~isempty(obj.ridge.gamma)
        fprintf('--> Best JK score: %g\n', obj.ridge.best_score);
    end
    fprintf('--> GOF:\n');
    disp(gof);
    fprintf('\n');
end


%
if ~isempty(fignum)
    t = linspace(0, 1e-3*binwidth*size(X_test0,2), size(X_test0,2));    
    nolabels = true;
    if train_drr == test_drr
        ysub = 2;
    else
        ysub = 3;
    end
    
    % Plot #1    
    figure(3 + fignum);
    clf;
    subplot(ysub,1,1);
    ax = spec.plot_spectrogram(t, 1e-3*obj.f, X_test0,...
        'fignum', 3+fignum, 'nolabels', nolabels);
    set(ax(1), 'XTickLabel', '');
    ylabel( aux.ctitle(sprintf('$X_{train}$ (DRR: %d)', train_drr), '\\') );
    title_str1 = sprintf('Sort type: %s, Speaker %d, CC: %.3f, NMSE: %.3f (BW: %g)',...
        sort_type, test_speaker, gof.CC, gof.nmse, binwidth);
    title(title_str1);
    
    subplot(ysub,1,2);
    ax(2) = spec.plot_spectrogram(t, 1e-3*obj.f, obj.X_est, 3+fignum, nolabels);
    set(ax(2), 'XTickLabel', '');
    ylabel( aux.ctitle(sprintf('$\\hat{X}_{test}$ (DRR: %d)', test_drr), 'Frequency (kHz)' ));    

    if train_drr ~= test_drr
        xlabel(ax(2), '');
        subplot(ysub,1,3);
        ax(3) = spec.plot_spectrogram(t, 1e-3*obj.f, X_test, 3+fignum, nolabels);
        ylabel( aux.ctitle(sprintf('$X_{test}$ (DRR: %d)', test_drr), '\\') );
    end
    xlabel('Time (sec)');
    linkaxes(ax);
    aux.abc(ax);
    
    % Plot #2
    figure(5 + fignum);
    clf;    
    n_lags = obj.n_lags;
    ax = obj.plot_reconstruction_filters('fontsize', 18, 'fignum', 5 + fignum);  
    drawnow;
end  


%% Plot arbitrary R-Filter
% Loads a struct with fields:
switch upper(data_type)
    case 'SU'
        units2plot =  [1 6 3 12 17 15];  % [2 6 11 15 10 20];    % SU
        
    case 'MUA'
        units2plot = [1 6 3 12 17 15];  %[1 14 3 12 17 15];  % MUA
end
        
xysub = [3 2];
obj.plot_reconstruction_filters('fontsize', 18, 'fignum', 199, 'units', units2plot, 'xysub', xysub);




