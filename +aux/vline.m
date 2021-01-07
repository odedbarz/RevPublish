function vline_h = vline(xpos, varargin)
% 
% function vline(plth, xpos)
% 
% Input:
%   ax:     (1x1) a handle to an axes to plot in;
%   xpos:   (1xN) position(s) to plot the vertical lines
%
% Description:
%   Plots vertical lines in ax.
% 

%% Parse the input
p = inputParser;

addRequired(p,'xpos', @(x) isnumeric(x) || iscell(x));

addOptional(p,'ax', [], @(x) ishandle(x) || isnumeric(x));
addOptional(p,'ZData', [], @isnumeric);

defaultColor = 'k';
addOptional(p,'color', defaultColor);

displayName = '';
addOptional(p,'displayName', displayName);

defaultLinestyle = '--';
addOptional(p,'LineStyle', defaultLinestyle);

defaultLineWidth = 2;
addOptional(p,'LineWidth', defaultLineWidth);

parse(p, xpos, varargin{:});
pars = p.Results;

xpos        = p.Results.xpos;
ax          = p.Results.ax;
linecolor   = p.Results.color;
linestyle   = p.Results.LineStyle;
displayName = p.Results.displayName;
lineWidth   = p.Results.LineWidth;

%%
if isempty(ax)
    ax = gca;
end

axes(ax);

dummy = ylim(ax);
y_min = dummy(1);
y_max = dummy(2);

hold on
% % Add the vertical lines
% vlineh = arrayfun(@(X) plot(X*[1, 1], [y_min, y_max],...
%     'Color', linecolor,...
%     'LineStyle', linestyle,...
%     'LineWidth', lineWidth,...
%     'DisplayName', displayName ),...
%     xpos, 'UniformOutput', false);

% Add the vertical lines
vline_h = [];
for q = 1:length(xpos)
    if iscell(xpos) 
        xq = xpos{q};
    elseif isvector(xpos)
        xq = xpos(q);
    else % its a scalar
        xq = xpos;
    end
    
    if isempty(pars.ZData)
        plth_q = plot(xq*[1, 1], [y_min, y_max]);
        vline_h = [vline_h, plth_q];
    else
        plth_q = plot3(xq*[1, 1], [y_min, y_max], 2*max(pars.ZData(:))*[1, 1]);
        vline_h = [vline_h, plth_q];
    end
end
hold off

ylim([y_min, y_max]);

set(vline_h,...
    'Color', linecolor,...
    'LineStyle', linestyle,...
    'LineWidth', lineWidth,...
    'DisplayName', displayName );




