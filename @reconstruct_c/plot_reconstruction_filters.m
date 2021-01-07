function ax = plot_reconstruction_filters(obj, varargin)
%
%   function ax = plot_reconstruction_filters(obj, ...)
%

p = inputParser;

% *** Required parameters ***
%addRequired(p, '', @isnumeric);      

% *** Optional parameters ***
addOptional(p, 'units', [], @isnumeric);                % filters to plot
addOptional(p, 'fignum', [], @isnumeric);               % 
addOptional(p, 'xysub', [], @isnumeric);                % 
addOptional(p, 'f', 1e-3*obj.f, @isnumeric);            % (Hz)
addOptional(p, 'fontsize', 18, @isnumeric);                % filters to plot
addOptional(p, 'max_subs_in_a_row', 5, @isnumeric);  	% 
addOptional(p, 'colormap_type', 'jet', @isstr);         % 
addOptional(p, 'add_colorbar', 0, @(x) islogical(x) || isscalar(x));        
addOptional(p, 'precision', 0.1, @isscalar);        

parse(p, varargin{:});
pars = p.Results;

assert(2 == length(pars.xysub) || 0 == length(pars.xysub),...
    '--> You need to provide [# of X, # of Y] subplot or just to leave it empty!');


%% Set the input
if isempty(obj.G)
    warning('--> There isn''t any filter(s) to plot! Use FIT to learn a new  one first!');
    return;
end

if ~isempty(pars.fignum)
    figure(pars.fignum);
    clf;
else
    pars.fignum = 99;
    figure(pars.fignum);
    clf;    
end


    
%%
% # of subplots
if isempty(pars.xysub) && 1 == length(pars.units)
    xsub = 1;
    ysub = 1;
elseif isempty(pars.xysub)
    xsub = min(obj.n_neurons, pars.max_subs_in_a_row);
    ysub = min( ceil(obj.n_neurons/xsub), pars.max_subs_in_a_row);
else
    xsub = pars.xysub(1);
    ysub = pars.xysub(2);
end

% Filters of units to plot
if isempty(pars.units)
    units = 1:(xsub*ysub);
else
    units = pars.units;
    units = units(1:min(length(units), xsub*ysub));
end
% n_neurons = min(obj.n_neurons, length(units));
n_neurons = length(units);

assert(max(units) <= obj.n_neurons, '--> You asked to plot UNITS that are NOT available!');


%% Plot
% the x-axis
x_lags = obj.lags * obj.binwidth;
% log_tick_fun = @(x) pars.precision*fix(1/pars.precision*10.^x);
ax = nan(1, xsub*ysub);
counter = 0;
c_axis = zeros(n_neurons, 2);

for k = 1:n_neurons
    if counter > xsub*ysub
        break;
    end
    
    counter = counter +1;
    if counter > xsub*ysub, break; end
    
    ax(k) = subplot(ysub, xsub, counter);
    
    % Get the k'th filter indices
    neuron_idx = units(k);
    filter_idx = obj.n_lags*(neuron_idx-1) + (1:obj.n_lags);        
    Gk = obj.G(filter_idx, :)';     % get the k'th filter
    
    if ~isvector(Gk)
        surf_h = surf(x_lags, pars.f, Gk);
        surf_h.LineStyle = 'none';
        %surf_h = contour(x_lags, log10(f), Gk);
        view(2);
    else
        % No lags; each reconstruction filter is a simple 1D FIR 
        surf_h = [];
        plot(1e-3*obj.f, Gk);
    end
    axis tight    
    %drawnow;
    
    %{
    yticks = get(ax(k), 'YTick');
    if 2<=length(yticks) 
        yticks = yticks(2:end-1);
        yticklabels = log_tick_fun(unique( [log10(pars.f(1)), yticks, log10(pars.f(end))] ));
    end
    set(ax(k), 'YTicklabels', yticklabels);
    %}
    set(ax(k), 'YDir', 'normal');
    colormap(pars.colormap_type);
    c_axis(k,:) = caxis;
    
    [Ik, Jk] = ind2sub([xsub,ysub], counter);
    
    if ~isempty(surf_h)
        aux.vline([0], 'ZData', surf_h.ZData);
    end
    
    % Set the y-axis
    if 1 == Ik && ysub == Jk
        ylabel(ax(k), 'Frequency (kHz)');     % y-axis    
        xlabel(ax(k), 'Lags (ms)');           % x-axis
    else
        set(gca, 'YTickLabel', '');
        set(gca, 'XTickLabel', '');
    end
end

%
ax = reshape(ax, [xsub, ysub]);
linkaxes(ax(1:n_neurons));
aux.abc(ax(1:n_neurons), 'label_type', 'numbers', 'fontsize', 45);

% Font size
for nn = 1:numel(ax)
    if isnan(ax(nn))
        continue;
    end
    set(ax(nn), 'FontSize', pars.fontsize);
end


% Add a colorbar and make sure that the axes doesn't change in size (push 
% the colorbar "outside").
if pars.add_colorbar
    % Set the CAXIS to be the same over ALL subplots
    caxis_new = [min(c_axis(:,1)), max(c_axis(:,2))];
    arrayfun(@(n) caxis(ax(n), caxis_new), 1:n_neurons );

    ax_j = ax(n_neurons);
    pos = get(ax_j, 'Position');
    colorbar(ax_j);
    set(ax_j, 'Position', pos);
end

% Add ABCD labels
%aux.abc(ax, 48, 'northwest');



