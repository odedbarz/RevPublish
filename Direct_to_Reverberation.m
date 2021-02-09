function [drr_values, direct_ending_sample] = Direct_to_Reverberation(h_rir, direct_ending_sample, fignum)
%
%   function [drr_values, direct_ending_sample] = Direct_to_Reverberation(h_rir, direct_ending_sample, fignum)
%
%
% Input:
%   h_rir: (Nx1) the overall reverberation impulse response(s)
%   [direct_ending_sample]: (optional 1x1) the ending (starts at the beginning)
%                       of the temporal window that defines the direct 
%                       response (default: 250 samples).
%
% Description: 
%   Calculates the Direct-to-reverberation ratio
%

if 2 > nargin || isempty(direct_ending_sample)
    direct_ending_sample = min(size(h_rir,1), 100);    % samples
end

if 3 > nargin
    fignum = [];
end

assert(1 == prod(max(h_rir)), '--> ERROR at [DRR.m]: all RIR columns must be normalize to 1!!');

% The direct impulse
h_direct = h_rir(1:direct_ending_sample,:);

% The reverberation part
h_reverb = h_rir((1+direct_ending_sample):end,:);

% Calc the DRR:
drr_values = nan(1, size(h_reverb,2));
for kk = 1:size(h_reverb,2)
    if sum(abs(h_reverb(:,kk))) > 1e-6
        direct_energy = 10*log10( sum(abs(h_direct(:,kk).^2)) );
        reverb_energy = 10*log10( sum(abs(h_reverb(:,kk).^2)) );
        drr_values(:,kk) = direct_energy - reverb_energy;

    else
        drr_values(:,kk) = Inf;
    end
end


%% DEBUG
if ~isempty(fignum)
    figure(fignum);
    clf;
    
    % ---------------
    ax = subplot(3,1,1);
    plth = plot(1:direct_ending_sample, h_direct(:,1), '--');
    hold on
    plth(end+1) = plot((1+direct_ending_sample):length(h_rir), h_reverb(:,1), 'Color', plth(1).Color);
    plot(direct_ending_sample*[1,1], ylim, '--k');    
    hold off       
    set(ax(1), 'XTickLabel', '');
    %xlabel('Samples');    
    legend(plth, 'Direct Response (Left)', 'Reverberation (Left)');
    xlim([1,10e3]);
    title(sprintf('Left Side, D/R: %s dB', num2str(drr_values(1), '%.2f  ')));

    % ---------------
    ax(2) = subplot(3,1,2);
    plth = plot(1:direct_ending_sample, h_direct(:,2), '--', 'Color', aux.rpalette('new02'));
    hold on
    plth(end+1) = plot((1+direct_ending_sample):length(h_rir), h_reverb(:,2), 'Color', plth(1).Color);
    plot(direct_ending_sample*[1,1], ylim, '--k');
    hold off       
    set(ax(2), 'XTickLabel', '');
    %xlabel('Samples');    
    legend(plth, 'Direct Response (Right)', 'Reverberation (Right)');
    title(sprintf('Right Side, D/R: %s dB', num2str(drr_values(2), '%.2f  ')));

    linkaxes(ax(1:2), 'x');
    xlim([1,10e3]);
    ylim([-0.1, 1.0]);
    
    % =============================
    %figure(10+fignum);
    %clf;    
    ax(3) = subplot(3,1,3);
    plot(1:length(h_rir), 20*log10(max(eps,h_rir)), '.');
    xlabel('Samples');    
    ylabel('$20log_{10}$');
    aux.hline(-60);
    x_text = 0.9*max(xlim);
    txth = text(x_text, -120, '$T_{60}$');
    set(txth, 'FontSize', 36)
end







