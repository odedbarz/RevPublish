function [ax_sub] = plot_FRA(S, varargin)
%
%   function plot_FRA(S, [figh], [syncchan], [spikechan])
%
% Input:
%   S         : (struct) An Impale's data structure of the FRA measurement.
%   [figh]    : (1x1) graphic handle)
%   [syncchan]: (1x1) the syncchan, usually equals 0.
%
% Description:
%   Plots the FRA(s), a line around the prominent active area(s), and the CF(s)
%   of that(these) area(s); the CF(s) is(are) defined as the minimum point(s) 
%   in this(these) area(s).
%

%% Set the inputs
p = inputParser;

addRequired(p, 'S', @(x) isstruct(x));           

addOptional(p, 'figh', [], @ishandle); 
% addOptional(p, 'ax', [], @ishandle); 
addOptional(p, 'syncchan', 0, @isnumeric);     
addOptional(p, 'spikechan', -1, @isnumeric);     

parse(p, S, varargin{:});

syncchan = p.Results.syncchan;
spikechan = p.Results.spikechan;
if isempty(p.Results.figh), figh = figure(99); else, figh = p.Results.figh; end



%%
fonts.title   = 16;
fonts.axes    = 14;
fonts.text    = 18;
xtick_precision = 2;

figure(figh);
clf;

% Get all recorded spikechan 
if -1 == spikechan
    spikechan = spiketools.get_all_spikechan(S, syncchan);
    if iscell(spikechan), spikechan = spikechan{:}; end
end

% Calc FRA for each of the available spikes
fra = arrayfun(@(SPK) medit.FRA(S, SPK, syncchan), spikechan, 'UniformOutput', false) ;

% Set the subplot dimensions
len_fra = length(fra);
sub_ratio = min(3, len_fra);
subx = ceil(len_fra/sub_ratio);
suby = sub_ratio;
ax_sub = nan(1, len_fra);


for kk = 1:len_fra
    ax_sub(kk) = subplot(subx,suby,kk);
    
    % Frequency
    innerSeq.f0  = S.innerSeq.master.values;   
    innerSeq.var = S.innerSeq.master.var;       
    
    % dB SPL
    outerSeq.db  = S.outerSeq.master.values;   
    outerSeq.var = S.outerSeq.master.var;       
    
    x = 1e-3*(innerSeq.f0+eps);
    y = outerSeq.db;
    z = fra{kk}.rates';
    
    % Set the image:
    surf_h = surf(ax_sub(kk), x, y, z);      % [kHz x SPL x FRA]; 

    view(ax_sub(kk), 2);
    set(gca, 'XScale', 'log');
    set(surf_h, 'LineStyle', 'none');
    colormap hot;
    axis tight;
    h = colorbar;
    h.FontSize = fonts.axes;

    xtick = round(10^xtick_precision*exp([log(0.1):log(18)]))/10^xtick_precision;   % [kHz]
    set(gca, 'XTick', xtick);

    spikechan_kk = spikechan(kk);
    channel_number_kk = viewer.spikenum_to_Enum(spikechan_kk);
    title_kk = sprintf('FRA E$_{%d}$ (CF: %g Hz, spikechan: %d)', channel_number_kk, fra{kk}.CF, spikechan_kk);
    title_h = title(title_kk, 'Interpreter', 'latex');
    set(title_h, 'Fontsize', fonts.title, 'FontWeight', 'bold');

    xlabel(ax_sub(kk), innerSeq.var, 'Interpreter', 'latex');
    ylabel(ax_sub(kk), outerSeq.var, 'Interpreter', 'latex');
    
    % Add the CF (kHz)
    xc = fra{kk}.graphics.xc;
    yc = fra{kk}.graphics.yc;
    zc = fra{kk}.graphics.zc;
    
    hold on
    plot3(xc, yc, zc, '--w', 'linewidth', 3);
    hold off
    
    if isempty(xc) || isempty(yc)
        continue;
    end
    
    % Add the CF TEXT
    CF = fra{kk}.CF;
    min_yc_idx = fra{kk}.graphics.min_yc_idx;    
    
    % Make sure the text is INSIDE the axes
    % Extent: [left bottom width height]    
    %xtext = 0.5*xc(min_yc_idx);
    %ytext = 0.95*yc(min_yc_idx);
    xtext = 0.8*exp(mean(log(xlim)));
    ytext = mean(ylim);
    
    ztext = zc(1);
    cf_text = [num2str(CF), ' Hz'];
    txt_h = text(xtext, ytext, ztext, cf_text, 'color', 'w');
    set(txt_h, 'FontSize', 22);
    xlabel('$F_{tone}$ [kHz]');
    ylabel('Amp. [dB SPL]');
    
    
    % Add a marker at the CF (minimum) point
    hold on
    plot3( xc(min_yc_idx), yc(min_yc_idx), zc(min_yc_idx), 'sw', 'linewidth', 4, 'markersize', 14);
    hold off
end




