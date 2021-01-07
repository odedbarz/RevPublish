function rasters_h = plot_dotraster(y_bias, t, ch, varargin)
%
%    function rasters_h = tps_trials = plot_dotraster(t, ch, y_bias, spikechan, syncchan, [dots_color])
%


%% Set the input
p = inputParser;

addRequired(p, 'y_bias', @isnumeric);               
addRequired(p, 't', @isnumeric);               
addRequired(p, 'ch', @isnumeric);               

addOptional(p, 'spikechan', 1, @isnumeric);        
addOptional(p, 'syncchan', 0, @isnumeric);        
addOptional(p, 'dots_color', [], @(x) isnumeric || ischar(x));        
addOptional(p, 'markersize', 6, @(x) isnumeric );        
addOptional(p, 'raster_type', 'bars', @ischar );        % {'dots', 'lines', 'bars'}     
addOptional(p, 'time_units', 'ms', @ischar );           % {'ms', 'sec'}
addOptional(p, 'duration_sec', [], @isnumeric );
addOptional(p, 'trial_to_plot', [], @isnumeric );       % (1x1) trial number to plot 

parse(p, y_bias, t, ch, varargin{:});

spikechan  = p.Results.spikechan;              
syncchan   = p.Results.syncchan;              
dots_color = p.Results.dots_color;              
markersize = p.Results.markersize;              
raster_type= p.Results.raster_type;              
time_units = p.Results.time_units;              
duration_sec= p.Results.duration_sec; 
trial_to_plot= p.Results.trial_to_plot; 

if isempty(t) 
    rasters_h = [];
    return;
end

% assert(fix(y_bias) == y_bias);      % y_bias must be an integer 
assert(length(t) == length(ch));

if strcmpi('sec', time_units)
    tu = 1e-3;  % ms to sec
else
    tu = 1;  % ms; don't change the units
end


%% Set the spike trains
% Get spike times for Spch stimulus
[tps, ch_new, ~] = spiketools.get_spike_times(t, ch, spikechan, syncchan);

% Divide the TPS (all recorded spikes) into trials
n_trials   = nnz(syncchan == ch_new);
assert(0 < n_trials);
trials     = cumsum( ch_new == syncchan );
tps_trials = arrayfun(@(IDX) tps(IDX==trials), 1:n_trials, 'UniformOutput', false );

% Used only with raster_type == 'bars'
if isempty(duration_sec)
    nbins = ceil(1.2*max(tps));             % samples
else
    nbins = 1:0.5:ceil(1e3*duration_sec);   % (ms) dt= 0.5
end

% Remove all-empty trials
is_empty_trial = cellfun(@(X) isscalar(X) && (0 == X), tps_trials, 'UniformOutput', true);
tps_trials = tps_trials(~is_empty_trial);

% Update the number of syncs
n_trials = length(tps_trials);


%% Plot the dot-raster
% All trials (N_SYNC) are plotted in on the interval Y_BIAS + [0,1] 
y_add = (1-0)/n_trials;
% rasters_h = zeros(1, n_trials);
rasters_h = cell(1, n_trials);

if isnan(y_bias)
    % Add no bias in this case
    y_bias_all = 0;
else
    y_bias_all = y_bias -0.5 + 0.5*y_add;
end

if isempty(trial_to_plot)
    trial_to_plot = 1:n_trials;
end

% hold on
bias_counter = 0;
for k = trial_to_plot
    bias_counter = bias_counter + 1;

    switch lower(raster_type)
        case 'dots'
            % dots
            x = tu * tps_trials{k};
            y = y_bias_all + y_add*(bias_counter-1) + zeros(1, length(x));
            rasters_h{k} = plot(x, y, '.'); 

        case 'lines'
            % Option #1: this option is slow & doesn't work so good!
            % Get color automatically from MATLAB
            x = tu * tps_trials{k};
            
            y_bias_lines = 0.0; % -0.15;
            y = y_bias_lines + y_bias_all + y_add*(bias_counter-1) + zeros(1, length(x));
            
            h = plot(x(1), y(1), '.'); 
            dots_color = h.Color;
            delete(h);

            % lines
            line_height = 0.1;
            x_ = [x(:)'; x(:)'];
            y_ = [y(:)' - line_height*y_add; y(:)' + line_height*y_add];
            rasters_h{k} = plot(x_, y_, '-', 'Color', dots_color); 

        case 'bars'
            [psth, bins] = hist(tps_trials{k}, nbins);
            rasters_h{k} = plot(tu * bins, y_bias_all + y_add*(bias_counter-1) + psth*0.6*y_add);
            set(rasters_h{k}, 'LineWidth', 0.1);

    otherwise
        error('Invalid RASTER_TYPE!!');
    end

    if ~isempty(dots_color)
        set(rasters_h{k}, 'Color', dots_color);
    else
        dots_color = get(rasters_h{k}, 'Color');
    end
    set(rasters_h{k}, 'markersize', markersize);
    
    hold on
end
hold off



