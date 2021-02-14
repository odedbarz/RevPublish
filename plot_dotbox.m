function pars = plot_dotbox(X, varargin)
%
% pars = plot_dotbox(X, [x_axis_jitter_std], [labels], [median_line_length], [median_line_width])
%

p = inputParser;
addRequired(p, 'X', @isnumeric);	% (n_neurons x n_cases)       

addOptional(p, 'x_axis_jitter_std', 0.05, @isnumeric);      % (1x1) jitter along the x-axis      
addOptional(p, 'labels', {}, @(x) iscell(x) || isnumeric(x));                      % (1 x n_cases)           
addOptional(p, 'median_line_length',  0.30, @isnumeric);    % (1x1)           
addOptional(p, 'median_line_width', 3, @isnumeric);         % (1x1)      
addOptional(p, 'marker_size',  6, @isnumeric);             % (1x1)      
addOptional(p, 'std_outliers_factor',  [], @isnumeric);   	% (1x1)      
addOptional(p, 'color', [], @(x) isnumeric(x) || isstr(x));   	% (1x1) adds a boxplot

% addOptional(p, 'verbose', 0, @isnumeric);           
addOptional(p, 'fignum', [], @isnumeric);     
parse(p, X, varargin{:});
pars = p.Results;


%%
% x_axis_jitter_std  = p.Results.x_axis_jitter_std;
labels             = p.Results.labels;
% median_line_length = p.Results.median_line_length;
% median_line_width  = p.Results.median_line_width;
% marker_size        = p.Results.marker_size;
% std_outliers_factor= p.Results.std_outliers_factor;
% verbose           = p.Results.verbose;
% fignum             = p.Results.fignum;

if ~isempty(pars.fignum)
    figure(pars.fignum);
    clf;
end

[n_neurons, n_cases] = size(X);
assert( (isempty(labels)) || n_cases == length(labels));



%% Stats
pars.median.fun = @(x) nanmedian(x);
pars.median.M = nan(1,n_cases);

pars.se.fun = @(x) mad(x);
pars.se.M = nan(1,n_cases);

% Outliers
if ~isempty(pars.std_outliers_factor)
    % Remove outliers
    I = (X - median(X,1)) <= pars.std_outliers_factor*mad(X(:),1);
else
    % Use all points; don't account for outliers
    I = true(size(X));
end



%% Plot the jittered dots
pars.ax = gca;
pars.points_h = zeros(1, n_cases);
pars.outliers_h = zeros(1, n_cases);
pars.outliers_markersize = 5*pars.marker_size;
pars.outliers_alpha = 0.7;
pars.outliers_linewidth = 1.5;
pars.median_h = zeros(1, n_cases);

max_number_of_colors = 7;

for k = 1:n_cases    
    if 1 == n_neurons
        x_axis_jutter_k = 0;
    else
        x_axis_jutter_k = pars.x_axis_jitter_std * randn(1,n_neurons);
    end
    
    % Colors:
    if isempty(pars.color)
        color_k = aux.rpalette(sprintf('new%02d', 1+mod(k-1,max_number_of_colors)));
    else
       color_k = pars.color; 
    end
    
    
    % Add a bit of a jitter
    x_k = k + x_axis_jutter_k;    
    
    % Plot the dots
    pars.points_h(k) = plot(x_k(I(:,k)), X(I(:,k),k),...
        'o',...
        'MarkerEdgeColor', color_k,...
        'MarkerFaceColor', color_k );
    set(pars.points_h(k), 'MarkerSize', pars.marker_size);            
    
    % Add the k'th median  
    pars.median.M(k) = pars.median.fun( X(:,k) );
    pars.se.M(k) = pars.se.fun( X(:,k) );
    hold on
    pars.median_h(k) = plot(k+pars.median_line_length*[-1, 1], pars.median.M(k)*[1, 1], 'Color', color_k);
    set(pars.median_h(k), 'LineWidth', pars.median_line_width);
    
    % Add the outliers (samples that do not contribut to the median
    if 0 < nnz(~I(:,k))
        pars.outliers_h(k) = scatter(x_k(~I(:,k)), X(~I(:,k),k), ...
            pars.outliers_markersize, ...
            'Marker', 'x',...    
            'MarkerEdgeColor', get(pars.points_h(k), 'Color'),...
            'MarkerEdgeAlpha', pars.outliers_alpha,...
            'LineWidth', pars.outliers_linewidth );
    end
    
end
hold off
axis tight


% Set the labels (if provided)
if ~isempty(labels)   
    if isnumeric(labels)
        % In this case convert the labels to a cell array
        labels = arrayfun(@(x) sprintf('%g',x), labels, 'uni', 0);
    end
    
    set(gca, 'XTick', 1:n_cases);
    set(gca, 'XTickLabel', labels);
end














