function [ax, surf_h] = plot_spectrogram(t, f, Sft, varargin) 
%
%   function surf_h = plot_spectrogram(t, f, Sft, [fignum], [nolabels], 
%       [fontsize], [log_scale]) 
%
% t  : (sec) time axis (n_smp)
% f  : (Hz) frequency axis (n_bands)
% Sft: (n_bands x n_smp) spectrogram to plot 
%
% surf_h (1x1): handle to the 
%
% Description:
%   Plots the spectrogram Sft. The function uses SURF (instead of IMAGSC), so
%   that log-scale frequencies are represented reliably along the y-axis.
%

p = inputParser;
addRequired(p, 't', @isnumeric);             % (ms)        
addRequired(p, 'f', @isnumeric);             % (kHz)        
addRequired(p, 'Sft', @isnumeric);           % (fxt) spectrogram       

addOptional(p, 'fignum', [], @(x) isnumeric(x) || ishandle(x)); 
addOptional(p, 'nolabels', 0, @(x) isnumeric(x) || islogical(x)); 
addOptional(p, 'fontsize', 14, @isnumeric); 
addOptional(p, 'log_scale', 1, @isnumeric); 
addOptional(p, 'precision', 4, @isnumeric);     % 10^precision

parse(p, t, f, Sft, varargin{:});
pars = p.Results;

if isempty(pars.fignum) || pars.fignum == 0
    % Create a new & clean figure
    pars.fignum = figure;
    figure(pars.fignum);
    clf;
end


%% Plot the spectrogram on log-scaled axes (using surf)
figure(pars.fignum);
ax = gca;

if pars.log_scale
    surf_h = surf(t, log2(f), Sft, 'LineStyle', 'none');
    view(2); 
    axis tight;
    set(ax, 'FontSize', pars.fontsize);
    yticks = get(ax, 'YTick');
    set(ax, 'YTickLabel', arrayfun(@(x) num2str(2.^x, pars.precision), yticks, 'uni', 0) );
else
    surf_h = surf(t, f, Sft, 'LineStyle', 'none');
    view(2); 
    axis tight;
    set(ax, 'FontSize', pars.fontsize);
end
drawnow;

if ~pars.nolabels
    xlabel('Time (sec)');
    ylabel('Frequency (kHz)');
    colorbar;
end

% surf_h.LineStyle = 'none';
set(ax, 'YDir', 'normal');
colormap jet


