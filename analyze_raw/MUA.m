function [y, pars] = MUA(x, fs, fs_dwn, fignum)
%
%   function [y, pars] = MUA(x, fs, fs_dwn, fignum)
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


%% Butterworth bandpass filter
bpf.order   = 3;
bpf.fcutlow = 300;  % (Hz)
bpf.fcuthigh= 4500; % (Hz)
[bpf.b, bpf.a] = butter(bpf.order, [bpf.fcutlow, bpf.fcuthigh]/(fs/2), 'bandpass');
bpf.x = filtfilt(bpf.b, bpf.a, x);


%% Clipping extreme values
% %{
clip.SD = std(bpf.x);
clip.thr= 2.0 * clip.SD;

switch pars.clipping_type
    case 0
        % (0) do nothing
        clip.x = bpf.x;

    case 1
        % (1) Hard clipping:
        clip.x = bpf.x .* ( abs(bpf.x) <= clip.thr );   % Stark et al., 2007

    case 2
        % (2) Soft-clipping
        clip.x = tanh(bpf.x./clip.thr).*clip.thr;

    otherwise
        error('--> [ERROR at MUA.m]: invalid pars.clipping_type!');

end
%}
% clip.x = bpf.x;   '!! NO CLIPPING !!'



%% Plot\DEBUG
if ~isempty(fignum)
    len_x = length(x);
    t = linspace(0, 1/fs*len_x, len_x)';     % (sec)
    f = linspace(0, fs, len_x);

    figure(fignum);
    fignum = fignum + 2;
    clf;
    ax = subplot(2,1,1);
    plot(t, [x(:,1), bpf.x, clip.x]);
    aux.hline(clip.thr(1)*[1, -1]);     % adds horizontal lines to the plot
    legend('$x$', '$x_{lpf}$', sprintf('$x_{clipped}$ (type: %d)', pars.clipping_type));
    xlabel('Time (sec)');
    ylabel('Amp.');
    %xlim([1.0, 1.25]);     % arbitrary section
    xlim([0.5954 0.7720]);  % arbitrary section

    ax(2) = subplot(2,1,2);
    plot(1e-3*f, db(abs(fft([x(:,1), bpf.x, clip.x]))));
    xlabel('Frequency (kHz)');
    ylabel('Mag. (dB)');
    xlim([0, 1e-3*fs/2])
end



%% RMS
nnl.x = clip.x.^2;  % full cutoff

% Lowpass
nnl.order = 2;
nnl.fc = 0.4*fs_dwn;       % (Hz) cuttoff frequency for the lowpass filter
[nnl.b, nnl.a] = butter(nnl.order, (nnl.fc/fs), 'low');
nnl.x_lpf = filtfilt(nnl.b, nnl.a, nnl.x);

% Downsample to 500 Hz
nnl.fs_new = pars.fs_dwn;
[P, Q] = rat(nnl.fs_new/fs);
nnl.xd = resample(nnl.x_lpf, P, Q);

% [MUA > 0] Make sure that there are no negative values
nnl.xd = max(0,nnl.xd);
y = sqrt(nnl.xd);

%     nnl.fs_new = fs;
%     y = nnl.x_lpf;

len_y = length(y);
t_dwn = linspace(0, 1/nnl.fs_new*len_y, len_y)';     % (sec)
f_dwn = linspace(0, nnl.fs_new, len_y)';

pars.f_dwn = f_dwn;
pars.t_dwn = t_dwn;


%% Plot\DEBUG
if ~isempty(fignum)
    len_x = length(x);
    t = linspace(0, 1/fs*len_x, len_x)';     % (sec)
    f = linspace(0, fs, len_x);

    zca = @(x) (x-mean(x))./std(x);

    figure(fignum);
    fignum = fignum + 2;
    clf;
    %ax = subplot(2,1,1);
    plot(t, zca(x(:,1)));
    hold on
    plot(t_dwn, zca(y(:,1)));
    hold off
    legend('$x$ (ZCA)', '$MUA$ (ZCA)');
    xlabel('Time (sec)');
    ylabel('Amp.');
    %xlim([1.0, 1.25]);     % arbitrary section

    figure(fignum);
    %fignum = fignum + 2;
    %ax(2) = subplot(2,1,2);
    plot(f_dwn, db(abs(fft(y))));
    xlabel('Frequency (Hz)');
    ylabel('Mag. (dB)');
    legend('$MUA$');
    xlim([0, pars.fs_dwn/2])

end
