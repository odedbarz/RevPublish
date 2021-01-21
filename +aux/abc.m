function [label_h, ax] = abc(ax, varargin)
%
%   function text_h = abc(ax, [fontsize], [location], [labels])
%
% 
% Valid supported Locations (as in legend):
%     * ['northeast']
%     * 'northeast'
%     * 'southeast'
%     * 'southwest'
%     * 'northwestoutside'
%     * 'northeastoutside'
%
% Notes:
% For a nice explanation, see
%   https://www.mathworks.com/matlabcentral/answers/286003-annotation-box-left-corner-position
% 
%

%% Parse the input
p = inputParser;

addRequired(p,'ax', @(x) all(ishandle(x(:))));

addOptional(p,'fontsize', 60, @isscalar);
addOptional(p,'location', 'northwest', @isstr);
addOptional(p,'labels', []);
addOptional(p,'label_type', 'uppercase', @isstr);
addOptional(p,'outside_xbias', 0.22, @isscalar);
addOptional(p,'outside_ybias', 0.05, @isscalar);

parse(p, ax, varargin{:});
pars = p.Results;
pars.location = lower(pars.location);

% label type
switch pars.label_type
    case 'uppercase'
        labels = arrayfun(@(x) x, 'abcdefghijklmnopqrstuvwxyz', 'uni', 0);  % lowercase
        labels = upper(labels);

    case 'lowercase'
        labels = arrayfun(@(x) x, 'abcdefghijklmnopqrstuvwxyz', 'uni', 0);
        
    case 'numbers'
        labels = arrayfun(@(x) sprintf('%d',x), 1:26, 'uni', 0);
        
    case 'romans'
         labels = {'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII'};
            
    case 'user'
        labels = pars.labels;
        
    otherwise
        error('--> ERROR @ [aux.abs.m]: unrecognized label type!!!');
end

assert(length(ax(:)) <= length(labels), '--> ERROR @ [aux.abs.m]: you have more axes than labels!')


%%
drawnow;    % make sure that all the rendering is over
figh = gcf;

ax_shape = size(ax);
ax   = ax(:);
n_ax = length(ax);

label_h = cell(1, n_ax);
label_counter = 0;  % in case of ~ishandle case

for k = 1:n_ax
    if ~ishandle(ax(k))
        continue;
    end

    set(ax(k), 'Units', 'normalized');
    
    label_counter = 1 + label_counter;
    %pos = get(ax(k), 'Position');
    pos = aux.plotboxpos(ax(k));
        
    % Dimension of the text box
    dim = pos.*[1 1 1 1];   % [left bottom width height]
        
    % Vertical alignment
    if contains(pars.location, 'north')
        vertical_alignment  = 'top';
    elseif contains(pars.location, 'south')
        vertical_alignment  = 'bottom';
    else
        error('--> Unrecognized location (%s)!', pars.location);
    end
    
    % Horizontal alignment
    if contains(pars.location, 'west')
        horizontal_alignment  = 'left';
        
        % Outside the axes?
        if contains(pars.location, 'outside')
            % Dimension of the text box
            dim(1) = dim(1) - pars.outside_xbias*dim(3);   % [left bottom width height]
            dim(1) = max(0, dim(1));
            dim(2) = dim(2) + pars.outside_ybias;      % good just for the "north" option...            
        end
        
    elseif contains(pars.location, 'east')
        horizontal_alignment  = 'right';
        
        % Outside the axes?
        if contains(pars.location, 'outside')
            % Dimension of the text box
            dim(1) = dim(1) + pars.outside_xbias*dim(3);   % [left bottom width height]
            dim(1) = min(1, dim(1));
            dim(2) = dim(2) + pars.outside_ybias;      % good just for the "north" option...            
        end
        
    else
        error('--> Unrecognized location (%s)!', pars.location);
    end
    
    % Set the text annotation
    label_str  = labels{label_counter};
    label_h{k} = annotation(figh,...
                            'textbox', dim,...
                            'String', label_str, ...
                            'FitBoxToText', 'On',...
                            'VerticalAlignment', vertical_alignment, ...
                            'HorizontalAlignment', horizontal_alignment,...
                            'FontSize', pars.fontsize, ...
                            'EdgeColor', 'none', ...
                            'Interpreter', 'latex' ...
                        );
end

ax      = reshape(ax, ax_shape);
label_h = reshape(label_h, ax_shape);

if 1 == n_ax
    label_h = label_h{1};
end




