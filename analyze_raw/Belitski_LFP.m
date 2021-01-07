function [y, Fs_dwn, t_dwn] = LFP(x, Fs, Fs_new, fignum)
%
%   function [y, Fs_dwn, t_dwn] = LFP(x, Fs, Fs_new, fignum)
%
% Description: Calculates the local field potential (LFP) in x.
%
% Reference:
% (1) Belitski, A., Gretton, A., Magri, C., Murayama, Y., Montemurro, M.A., 
%     Logothetis, N.K. and Panzeri, S., 2008. Low-frequency local field 
%     potentials and spikes in primary visual cortex convey independent 
%     visual information. Journal of Neuroscience, 28(22), pp.5696-5709.
%

% In Impale, Fs = 100k Hz by default
if 2 > nargin 
    Fs = 100e3;     % (Hz)
end

if 3 > nargin || isempty(Fs_new)
    Fs_new = 500;     % (Hz)
end

if 4 > nargin
    fignum = [];     %
end

% downsampling factors
[N, D] = rat(Fs/Fs_new);

% Remove the DC
x_ = x - mean(x);

% Apply a cos window on the boundaries to avoid discountiuous jumps when
% filtering
win.ms = 200;   % (ms)
win.smp = ceil(1e-3*win.ms * Fs);
win.Hon = cos(2*pi*linspace(-1/4,0,win.smp))';    % onset
win.Hoff = cos(2*pi*linspace(0,1/4,win.smp))';     % offset
x_(1:win.smp) = win.Hon.*x_(1:win.smp);
x_(end-[win.smp:-1:1]-1) = win.Hoff.*x_(end-[win.smp:-1:1]-1);


%% LPF 
d_lpf = designfilt('lowpassiir', 'FilterOrder', 8, ...
         'PassbandFrequency', 1e-3*Fs_new/2,... %0.225,...   % (kHz) Nyquist frequency
         'PassbandRipple',0.01, ...
         'StopbandAttenuation', 60, ...
         'SampleRate', 100);

x_lpf = d_lpf.filtfilt(x_);

if fignum
    figure(fignum);
    fignum = fignum + 3;
    clf;
    f = 1e-3*linspace(0, Fs, length(x_lpf))';
    plot(f, db(abs(fft([x_, x_lpf]))));
    legend('$x_{nnl}$', '$x_{lpf}$');
    xlabel('Frequency (kHz)');
    ylabel('dB');
end


%% Envelope
x_env = x_lpf;  % do none; keep all info (envelope + phase)
%x_env = abs(hilbert(x_lpf));

if fignum
    figure(fignum);
    fignum = fignum + 3;
    clf;
    dt = 1/Fs;              % (sec)
    t = linspace(0, dt*length(x_lpf), length(x_lpf));    
    plot(t, [x_lpf, x_env]);
    legend('$x_{lpf}$', '$x_{env}$');
    xlabel('Time (sec)');
    ylabel('Amp');
end


%% Downsample
Fs_dwn = Fs*(D/N);
assert(norm(Fs_dwn-Fs_new)<=5*eps);

% Downsample from Fs to Fs_new
y = resample(x_env, D, N);

% Remove the mean ("filter" the DC index)
%y = y - mean(y);

dt_dwn = 1/Fs_dwn;      % (sec)
t_dwn = linspace(0, dt_dwn*length(y), length(y));    


if fignum
    % f-axis for Fs
    f = 1e-3*linspace(0, Fs, length(x_lpf))';
    dt = 1/Fs;              % (sec)
    t = linspace(0, dt*length(x_lpf), length(x_lpf));

    % f-axis downsampled
    f_dwn = linspace(0, f(end)/N, length(x_lpf)/N)'; 
    
    figure(fignum);
    fignum = fignum + 5;
    clf;    
    subplot(1,2,1);
    plot(f, db(abs(fft(x_lpf))));
    xlabel('Frequency (kHz)');
    ylabel('dB');    
    subplot(1,2,2);
    plot(t, x_lpf);
    legend('$x_{bpf}$');
    xlabel('Time (sec)');
    ylabel('Amp');
    
    figure(fignum);
    clf;    
    subplot(1,2,1);
    plot(f_dwn, db(abs(fft(y))));
    xlabel('Frequency (kHz)');
    ylabel('dB');
    subplot(1,2,2);
    plot(t_dwn, y)
    legend('$y$');
    xlabel('Time (sec)');
    ylabel('Amp (dB)');
end

























