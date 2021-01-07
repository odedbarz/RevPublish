function ax = errorbar(x, y, dy, c, alpha)
% 
% function errorbar(x, y, dy, c, alpha)
% 
% Example:
%     aux.errorbar(psth.bins, psth.psth{idx0}, psth.err{idx0}, [0.1 0.2 0.9])
% 

if 4 > nargin, c = [.9 .9 .9]; end
if 5 > nargin, alpha = 0.6; end
if 6 > nargin, lw = 2; end


ax.fill = fill([x; flipud(x)], [y-dy; flipud(y+dy)], c,...
    'linestyle', 'none');
ax.fill.FaceAlpha = alpha;

hold on
ax.line = plot(x, y, 'Color', c, 'LineWidth', lw);
hold off


