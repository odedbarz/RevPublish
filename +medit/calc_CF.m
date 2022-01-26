function [CF, graphics] = calc_CF(S, fra_rates, Nq, n_countors, p_thr, n_min_bars, fignum)
%
% function [CF, graphics] = calc_CF(S, fra_rates, Nq, [fignum])
%
% Input:
%   S        : (struct) An Impale's data structure of the FRA measurement.
%   fra_rates: (NxM) the (dB x Frequency) rates.
%   Nq       : (1x1) the interpulation factor.
%   [fignum] : (1x1) optional, for plotting
%
% Output:
%   CF       : (1x1) the characteristic frequency of the neuron.
%   graphics : (struct) holds data on the boundary line of neuron's most prominent region.
%
% Description:
%   Calculates the FRA region and finds the CF as the minimum point of this
%   area.
% 
%
%

if 3 > nargin || isempty(Nq)
    Nq = 10;    % the interpulation factor
end

if 4 > nargin
    n_countors = 10;
end

if 5 > nargin
    p_thr = 5; % for prctile
end

if 6 > nargin
    n_min_bars = 4; % for prctile
end

if 7 > nargin
    fignum = [];    % the interpulation factor
end


assert(strcmpi(S.innerSeq.master.var, 'Frequency'), '--> Error at [calc_CF.m]: cannot find the ''Frequency'' variable!');
F0 = S.innerSeq.master.values;     % [Hz]

assert(strcmpi(S.outerSeq.master.var, 'dB SPL'), '--> Error at [calc_CF.m]: cannot find the ''dB SPL'' variable!');
spl = S.outerSeq.master.values;     % [dB SPL]



%%
thr_low = prctile(fra_rates(:), p_thr);
thr_high = prctile(fra_rates(:), 100-p_thr);
m = median(fra_rates(:)); 
if abs(m-thr_low) < abs(m-thr_high)
    % Exitatory neuron
    fra_rates( thr_high > fra_rates ) = 0.0;
else
    % Inhibitory neuron
    fra_rates( thr_low < fra_rates ) = 0.0;
end



%%
if sum(fra_rates) == 0
    CF = nan;
    graphics = [];
    return
end


%%
hist_range = [min(fra_rates(:)) : max(fra_rates(:))];
if nnz( histcounts( fra_rates(:), hist_range ) ) <= n_min_bars
%     CF = nan;
%     graphics = nan;
%     return;
    fra_rates = 0 * fra_rates;
end



%%
[Nf0, Nspl] = size(fra_rates);
% assert( strfind(), '--> This is not a FRA measurement' )
assert( Nf0 == length(F0), '--> [FRA2CF.m]: ERROR! Nf0 ~= length(F0)' );
assert( Nspl == length(spl), '--> [FRA2CF.m]: ERROR! Nf0 ~= length(F0)' );

F = griddedInterpolant( fra_rates );

% interpulation values:
xq = 1:1/Nq:(Nf0);
yq = 1:1/Nq:(Nspl);
F0_q  = interp1( 1:Nf0, F0, xq );
spl_q = interp1( 1:Nspl, spl, yq );

% Interpulated FRA:
fra_q = F( {xq, yq} );


%% Calc. the CF
se = @(xx) strel('disk', xx);
Iobr = fra_q';
Iobrd = imdilate(Iobr, se(0.5*Nq));
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);

% get the contours:
C = contourc(1e-3*(F0_q+1e-6), spl_q, Iobrcbr, n_countors);

% Get the longest contour
Lc = nan(1,size(C,2));
Lc_idx = nan(1,size(C,2));
nn = 1;
counter = 0;
while nn < size(C,2)
    counter = counter+1;
    Lc(counter) = C(2,nn);
    Lc_idx(counter) = nn;
    nn = nn + Lc(counter) +1;        
end
Lc = Lc(1:counter); 
[max_len, max_idx] = max(Lc);
idx2line = Lc_idx(max_idx);

xc = C(1,idx2line+(1:max_len));
yc = C(2,idx2line+(1:max_len));
zc = 2*max(fra_q(:))*ones(size(xc));

% Text of the CF on the image:
[~,min_yc_idx] = min(yc);
% cf = exp(xc(min_yc_idx));
CF = xc(min_yc_idx); 	% [kHz]
CF = fix(1e3*CF);       % [kHz] --> [Hz]

% Save the graphics for plotting later
graphics.min_yc_idx = min_yc_idx;
graphics.xc = xc;
graphics.yc = yc;
graphics.zc = zc;


%% Plot the CF
if ~isempty(fignum)
    fontsize.title   = 20;
    fontsize.axes    = 14;
    fontsize.text    = 18;
    linewidth.ctext  = 2;   % contour's text
    linewidth.cmarker= 3;

    % contour properties:
    co.line_style    = '--w';
    co.text_style    = 'w';
    co.marker_style  = 'sw'; 

    %precision        = 3;    % num2str precision
    xtick_precision  = 2;

    figure(fignum);
    set(gca, 'FontSize', fontsize.axes);

    % Set the image:
    surf_h = surf(1e-3*(F0_q+1e-6), spl_q, fra_q');      % [kHz x SPL x FRA]; 

    view(2);
    set(gca, 'XScale', 'log');
    set(surf_h, 'LineStyle', 'none');
    axis tight;
    h = colorbar;
    h.FontSize = fontsize.axes;
    colormap hot;

    hold on
    plot3(xc, yc, zc, co.line_style, 'linewidth', linewidth.cmarker);
    hold off

    xtick = round(10^xtick_precision*exp([log(0.1):log(18)]))/10^xtick_precision;   % [kHz]
    set(gca, 'XTick', xtick);

    xtext = 0.5*xc(min_yc_idx);
    ytext = 0.95*yc(min_yc_idx);
    ztext = zc(1);
    cf_text = [num2str(CF), ' Hz'];

    txt_h = text(xtext, ytext, ztext, cf_text, 'color', co.text_style);
    set(txt_h, 'FontSize', fontsize.text);
    xlabel('$F_{tone}$ [kHz]');
    ylabel('Amp. [dB SPL]');
    hold on
    plot3( xc(min_yc_idx), yc(min_yc_idx), zc(min_yc_idx), co.marker_style, 'linewidth', linewidth.ctext);
    hold off

    title_h = title('FRA');
    set(title_h, 'Fontsize', fontsize.title, 'FontWeight', 'bold');
end


