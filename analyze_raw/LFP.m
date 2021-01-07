function [y, pars] = LFP(x, fs, fs_dwn, fignum)
%
%
%
% Reference:
% Stark, E. and Abeles, M., 2007. Predicting movement from multiunit activity. 
% Journal of Neuroscience, 27(31), pp.8387-8394.
%

if 3 > nargin
    fs_dwn = 500;     % (Hz)
end

if 4 > nargin
    fignum = [];     % figure number
end

pars.fs_input = fs;         % (Hz) sampling rate of the input signal
pars.f_dwn = [];
pars.t_dwn = [];
pars.clipping_type = 2;
pars.fs_dwn = fs_dwn;       % (Hz) sampling rate of the output signal, after downsampling


%% RMS 
% nnl.x = clip.x.^2;      % full cutoff
nnl.x = x - mean(x);    % remove the mean (instead of a bandpass at 1Hz)
% nnl.x = x;              % do nothing


%% Butterworth bandpass filter
%{
bpf.order   = 2;
bpf.fcutlow = 1;    % (Hz) 
bpf.fcuthigh= 100;  % (Hz) Default: in Stark et al., it is 100 Hz
[bpf.b, bpf.a] = butter(bpf.order, [bpf.fcutlow, bpf.fcuthigh]/(fs/2), 'bandpass');
bpf.x = filtfilt(bpf.b, bpf.a, nnl.x);
%}

% Lowpass
bpf.order = 3;
bpf.fc = 200;       % (Hz) cuttoff frequency for the lowpass filter
[bpf.b, bpf.a] = butter(bpf.order, (bpf.fc/fs), 'low');
bpf.x = filtfilt(bpf.b, bpf.a, nnl.x);




%% Downsample to 500 Hz
dwnsmp.fs_new = pars.fs_dwn;
[P, Q] = rat(dwnsmp.fs_new/fs);
dwnsmp.xd = resample(bpf.x, P, Q);


%% Clear power-grid harmonics
f_ = linspace(0, dwnsmp.fs_new, size(dwnsmp.xd,1))';
% [~, idx60] = min(abs(f_ - 60));
[~, idx120] = min(abs(f_ - 2*60));
% [~, idx180] = min(abs(f_ - 3*60));
% [~, idx240] = min(abs(f_ - 4*60));

Xd = dct(dwnsmp.xd);
Xd(idx120 + [-1:1],:) = repmat(median(Xd(idx120 + [-3:3],:)), 3, 1);
xd_filtered = idct(Xd);

if ~isempty(fignum)     % DEBUG
    figure(fignum-1);
    clf;    
    plot(f_, db(abs(fft([dwnsmp.xd(:,1), xd_filtered(:,1)]))), '.-');
    xlim([0, 250]);
    aux.vline([60*[1, 2, 3, 4]], 'LineWidth', 0.5);
    xlabel('Frequency (Hz)');
    ylabel('Amp. (dB)');
    
    figure(fignum-2);
    clf;    
    plot([dwnsmp.xd(:,3) - xd_filtered(:,3)]);
    xlabel('Time (samples)');    
end


%%
% % [LFP > 0] Make sure that there are no negative values
% y = max(0, dwnsmp.xd);
% y = dwnsmp.xd;
y = xd_filtered;


len_y = length(y);
t_dwn = linspace(0, 1/dwnsmp.fs_new*len_y, len_y)';     % (sec)
f_dwn = linspace(0, dwnsmp.fs_new, len_y)';


pars.f_dwn = f_dwn;
pars.t_dwn = t_dwn;


%% Plot\DEBUG
if ~isempty(fignum)
    len_x = length(x);
    t = linspace(0, 1/fs*len_x, len_x)';     % (sec)
    f = linspace(0, fs, len_x);
    %t_dwn = linspace(0, 1/nnl.fs_new*len_x_dwn, len_x_dwn)';     % (sec)
    f_dwn = linspace(0, dwnsmp.fs_new, len_y)';
    
    figure(fignum);
    fignum = fignum + 2;
    clf;
    plot(1e-3*f, db(abs(fft([x, bpf.x]))));
    xlabel('Frequency (kHz)');
    ylabel('Mag. (dB)');
    legend('$x$', '$x_{lpf}$');
    xlim([0, 1e-3*pars.fs_input/2])
    
    figure(fignum);
    fignum = fignum + 2;
    clf;
    %ax = subplot(2,1,1);
    plot(t, x/std(x(:))*std(y));
    hold on
    plot(t_dwn, y);
    hold off
    legend('$x$ (normalized)', '$LFP$');
    xlabel('Time (sec)');
    ylabel('Amp.');    
    %xlim([1.0, 1.25]);     % arbitrary section
    
    figure(fignum);
    %fignum = fignum + 2;    
    %ax(2) = subplot(2,1,2);
    plot(f_dwn, db(abs(fft(y))));
    xlabel('Frequency (Hz)');
    ylabel('Mag. (dB)');
    legend('$LFP$');
    xlim([0, pars.fs_dwn/2])
    
end












