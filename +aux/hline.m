function hline_h = hline(ypos, varargin)
% 
% function hline(plth, ypos)
% 
% Input:
%   plth: (1x1) a handle of a plot object (plth = plot(...) )
%   ypos: (1xN) positions to plot the horizontal lines
%
% Description:
%   Plots vertical lines.
% 


%% Parse the input
p = inputParser;
addRequired(p,'ypos', @(x) isnumeric(x) || iscell(x));

addOptional(p,'ax', [], @(x) ishandle(x) || isnumeric(x));
addOptional(p,'ZData', [], @isnumeric);

defaultColor = 'k';
addOptional(p,'color', defaultColor);

displayName = '';
addOptional(p,'displayName', displayName);

defaultLinestyle = '--';
addOptional(p,'linestyle', defaultLinestyle);

defaultLineWidth = 2;
addOptional(p,'lineWidth', defaultLineWidth);

parse(p, ypos, varargin{:});

ypos        = p.Results.ypos;
ax          = p.Results.ax;
linecolor   = p.Results.color;
linestyle   = p.Results.linestyle;
displayName = p.Results.displayName;
lineWidth   = p.Results.lineWidth;
ZData       = p.Results.ZData;



%%
if isempty(ax)
    ax = gca;
end

axes(ax);

dummy = xlim(ax);
x_min = dummy(1);
x_max = dummy(2);

hold on
% % add vertical lines:
% hline_h = arrayfun(@(X) plot([x_min, x_max], X*[1, 1], 'Color', 'b'),...
%     ypos, 'UniformOutput', false);

% Add the vertical lines
hline_h = [];
for q = 1:length(ypos)
    if iscell(ypos) 
        xq = ypos{q};
    elseif isvector(ypos)
        xq = ypos(q);
    else % its a scalar
        xq = ypos;
    end
        
    if isempty(ZData)
        plth_q = plot([x_min, x_max], xq*[1, 1]);
        hline_h = [hline_h, plth_q];
    else
        plth_q = plot3([x_min, x_max], xq*[1, 1], 2*max(pars.ZData(:)));
        hline_h = [vline_h, plth_q];
    end
    
end
hold off


set(hline_h,...
    'Color', linecolor,...
    'LineStyle', linestyle,...
    'LineWidth', lineWidth,...
    'DisplayName', displayName );


% % Update all other fields of the lines:
% for kk = 1:length(ypos)    
%     hline_h{kk}.Color      = linecolor;
%     hline_h{kk}.LineStyle  = linestyle;
%     hline_h{kk}.DisplayName= displayName;
%     hline_h{kk}.LineWidth  = lineWidth;
% end
hold off


