function [ax, pars] = plot_strf(obj, add_margins, sigma)
%
% function [ax, surfh] = plot_strf(obj, [add_margins])
%
% 


%% Inputs
n = 1;  % every object holds ONE STFE only!
assert(isscalar(n), '--> [@strf_c]: n must be a scalar!');
if 1 >= nargin || isempty(add_margins)
    add_margins = false;
end
                
if 2 >= nargin
    sigma = [];
elseif 0 == sigma
    sigma = [];
end


%% Set the margins (1 or 3 axes)
if add_margins
    clf;
    ax(1) = subplot(10,10,[21,98]);     % MTF (2D)
    ax(2) = subplot(10,10,[1,18]);      % tMTF
    ax(3) = subplot(10,10,[29,100]);    % sMTF
else
    ax = gca;
end



%% *** Plot the STRF ***
% Time axis
tt = obj.lags * obj.binwidth;   % (ms)

axes(ax(1));
if isempty(sigma)
    strf = obj.strf;
else
    strf = imgaussfilt(obj.strf, sigma);
end

surfh = surf(tt, 1e-3*obj.f, strf );    % for the LOG scale!
set(ax, 'YDir', 'normal');
set(surfh, 'EdgeColor', 'none');
% shading interp
view(0,90);  % <--> view(0,90);     
ylabel('Frequency (kHz)');
xaxis_label = 'Lags (ms)';
xlabel(xaxis_label);


%% Add the BF mark to the plot
[BF, amp, bf_row_idx, bf_col_idx] = obj.calc_bf;

hold(ax(1), 'on');
plth = plot3(ax(1), tt(bf_col_idx), 1e-3*obj.f(bf_row_idx), 2*amp, 'xr',...
    'LineWidth', 6.0);


%% Add the temporal & spectral margins
if add_margins
    % *** Plot the FREQUENCY marginal (a function of time!) ***
    axes(ax(2));
    set(ax(2), 'XTickLabel', '');
    plot(tt, obj.strf(bf_row_idx,:));
    aux.hline(0);
    linkaxes(ax([1,2]), 'x');
    set(ax(2), 'XTickLabel', '');    

    % *** Plot the TEMPORAL marginal (a function of frequency!) ***
    axes(ax(3));
    
    % LEFT axes
    yyaxis('left');
    set(ax(3), 'YTickLabel', '');

    % RIGHT axes
    yyaxis('right');
    plot( obj.strf(:,bf_col_idx), 1e-3*obj.f, 'Color', aux.rpalette('new01'));
    aux.vline(0);    
    set(gca, 'YColor', [0,0,0]);
    axis(gca, 'tight');
    set(ax(3), 'YTickLabel', '');
    ax_title = 2;
    
else
    xlabel(ax, xaxis_label);
    ax_title = 1;
    colorbar(ax);
end

axis(ax, 'tight');
arrayfun(@(AX) hold(AX, 'off'), ax);
title(ax(ax_title), sprintf('$BF_{STRF}$: %.2f kHz', 1e-3*BF) );


pars.plth = plth;
pars.surfh= surfh;



