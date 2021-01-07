function hnd = twoSidedHorizontalArrow(xc, yc, width)

p1 = [xc-width/2, yc];  	% First Point
p2 = [xc+width/2, yc];  	% Second Point
dp = p2-p1;      % Difference

hold on
hnd(1) = quiver(p1(1), p1(2), dp(1), dp(2), 0);
hnd(2) = quiver(p2(1), p2(2), -dp(1), -dp(2), 0);
hold off

set(hnd, 'LineWidth', 2.0);
set(hnd, 'Color', 'k');
% set(hnd, 'AutoScale', 'on', 'AutoScaleFactor', 10);
set(hnd, 'MaxHeadSize', 0.5);

