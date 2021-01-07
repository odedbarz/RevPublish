% reconstruct_over_units.m
%
%
% Check reconstructions as a function of # of units

% Train: orthogonalize the measurements and get the most effective set of units
train_drr = 3     % drr.dry;
test_drr  = 3     % drr.dry;
n_units   = 10;    % Choose # of units to use for the reconstruction
% n_next    = 10;

drr = get_DRR_list_and_indices;

aux.cprintf('string', '*** SECTION #2 ***\n');
aux.cprintf('string', '--> train_drr   : %d\t (%s)\n', train_drr, drr.labels{train_drr});
aux.cprintf('string', '--> test_drr    : %d\t (%s)\n', test_drr, drr.labels{test_drr});

% Chhose the sorting algorithm:
% sort_type = 'ascend';
% sort_type = 'random';
sort_type = 'orthogonal';
% sort_type = 'invariant-strf';


H_train = squeeze( H_unit(:, train_drr, :) );

% Get an optimal set
switch sort_type
    case 'orthogonal'
        [slc.optimal_sorted, P, sv] = find_best_unit_set(H_train, 'n_svd', 36);
        
        % transpose to get the neurons for each orthogonal vector
        slc.optimal_sorted = unique(slc.optimal_sorted', 'stable');
        
    case 'random'
        slc.optimal_sorted = randperm(size(H_unit,3))';     % a random permutation
        %aux.cprintf('r', '--> optimal_sorted: RANDPERM\n');
        
    case 'ascend'
        slc.optimal_sorted = (1:size(H_unit,3))';     % a ascending order
        
    case 'invariant-strf'
        dummy = load('../.data/STRFs/invariant_bf_list.mat');
        slc.optimal_sorted = find(dummy.inv_bf);
        
        
    otherwise
        error('--> Unrecognized sort type!');
end
aux.cprintf('string', '--> sort_by    : %s\n', sort_type);

% Use these units
slc.optimal_sorted = slc.optimal_sorted(1:n_units, :);
slc.optimal_sorted = slc.optimal_sorted(:);     % make it a column vector
aux.cprintf('string', '--> # of units: %d\n', n_units);


% Training spectrogram
X1 = spec_st.Sft{train_drr};

% Testing spectrogram
X2 = spec_st.Sft{test_drr};

% units_list = 1:n_next:n_units;
units_list = 1:n_units;
sprintf('--> units_list: [%s]\n', [num2str(slc.optimal_sorted(1:n_units)', '%d ')])


clear gof
gof.CC  = nan(n_splits, n_units);
gof.mse = nan(n_splits, n_units);
gof.nmse= nan(n_splits, n_units);


for nn = 1:n_units
    n_units_used = units_list(nn);
    aux.vprint(verbose, '\n-> # of units: %d\n', n_units_used);

    % Use these units
    slc.optimal_units_used = slc.optimal_sorted(1:n_units_used);

    % Training responses
    y1 = squeeze( H_unit(:, train_drr, slc.optimal_units_used) );

    for sp = 1:n_splits
        aux.vprint(verbose && ~rem(sp-1, 3), '--> speaker: %d\n', sp);

        % *** Train ***         
        % Train the model using the "best" chosen units
        [X_train, X_test0, y_train, ~, ~] = train_test_split(X1, y1, ...
            'split_time_idx', split_time_idx,...
            'test_grp', sp );

        % Initialize the reconstruction object
        obj = reconstruct_c(binwidth,...
            'f', spec_st.f,...
            'iscausal', iscausal, ...
            'lags_ms', lags_ms ); 

        % Fit the model
        obj.fit(X_train, y_train,...
            'jk_flag', jk_flag,...
            'n_splits', n_splits, ...
            'fignum', []);


        % *** Test *** 
        % Testing responses
        y2 = squeeze( H_unit(:, test_drr, slc.optimal_units_used) );

        [~, X_test, ~, y_test, ~] = train_test_split(X2, y2, ...
            'split_time_idx', split_time_idx,...
            'test_grp', sp );

        % PREDICTION; Predict the spectrogram using the selected responses
        obj.predict(y_test);    

        % GOODNEES-of-FIT
        gof_n = goodness(X_test, obj.X_est);
        gof.CC(sp,nn)  = gof_n.CC;
        gof.mse(sp,nn) = gof_n.mse;
        gof.nmse(sp,nn)= gof_n.nmse;
    end
end



%% Plot stuff
% Plot the projections 
if ~isempty(fignum)
    %{
    figure(20+fignum);
    errorbar(1:n_splits, mean(gof.CC,2), std(gof.CC,[],2));
    xlabel('Speaker');
    xlim([0.5, 0.5+n_splits]);
    ylabel('CCs');
    title(sprintf('CCs vs. Speakers (averaged over units, BW: %g ms)', binwidth));
    %}
    
    figure(23+fignum);
    %clf;
    ax = errorbar(units_list, mean(gof.CC,1), std(gof.CC,[],1), ':.');
    set(ax, 'MarkerSize', 24, 'LineWidth', 2);    
    xlabel('Number of Units');
    xlim([0.5, 0.5+n_units]);
    ylim([0.3, 1.0])
    ylabel('CCs');
    title(sprintf('Sort Type: %s (averaged over speakers, BW: %g ms)', sort_type, binwidth));
    set(gca, 'FontSize', 24);
    legend_str = sprintf('train: %d (%s), test : %d (%s)', train_drr, drr.labels{train_drr},...
        test_drr, drr.labels{test_drr});
    legend(legend_str, 'Location', 'southeast');
    
    
    figure(25+fignum);
    clf;
    dot_pars = plot_dotbox(gof.CC, 'labels', units_list);
    ax = dot_pars.ax;
    xlabel('Number of Units');
    %xlim([0.5, 0.5+n_units]);
    %ylim([0.3, 1.0])
    ylabel('CCs');
    title(sprintf('Sort Type: %s (averaged over speakers, BW: %gms)', sort_type, binwidth));
    set(gca, 'FontSize', 24);
    legend_str = sprintf('train: %d (%s), test : %d (%s)', train_drr, drr.labels{train_drr},...
        test_drr, drr.labels{test_drr});
    legend(legend_str, 'Location', 'southeast');

    
    figure(27+fignum);
    clf;
    plot(units_list, gof.CC([1,4,7,12], :)', '.:');
    legend('speaker 1', 'speaker 4', 'speaker 7', 'speaker 12', 'Location', 'southeast');
    title(sprintf('Sort Type: %s (BW: %g ms)', sort_type, binwidth));
    xlabel('Number of Units');
    xlim([0.5, 0.5+n_units]);
    ylim([0.3, 1.0])
    ylabel('CCs');
    set(gca, 'FontSize', 24);
end



