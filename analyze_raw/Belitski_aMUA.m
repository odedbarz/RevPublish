function [y, Fs_dwn, t_dwn] = aMUA(x, Fs, Fs_new, fignum)
%
%   function [y, Fs_dwn, t_dwn] = aMUA(x, Fs, Fs_new, fignum)
%
% Description: Calculates the analog multi-unit activity (aMUA) in x.
%
% Reference:
% (1) Belitski, A., Gretton, A., Magri, C., Murayama, Y., Montemurro, M.A., 
%     Logothetis, N.K. and Panzeri, S., 2008. Low-frequency local field 
%     potentials and spikes in primary visual cortex convey independent 
%     visual information. Journal of Neuroscience, 28(22), pp.5696-5709.
% (2) Sadeghi, M., Zhai, X., Stevenson, I.H. and Escabí, M.A., 2019. A 
%     neural ensemble correlation code for sound category identification. 
%     PLoS biology, 17(10), p.e3000449.
% (3) Schnupp, J.W., Garcia-Lazaro, J.A. and Lesica, N.A., 2015. Periodotopy 
%     in the gerbil inferior colliculus: local clustering rather than a 
%     gradient map. Frontiers in neural circuits, 9, p.37.
%

% In Impale, Fs = 100k Hz by default
if 2 > nargin
    Fs = 100e3;     % (Hz)
end

if 3 > nargin || isempty(Fs_new)
    Fs_new = 1000;     % (Hz)
end

if 4 > nargin
    fignum = [];     %
end

filts.hpf = 'none';         % {'none', 'butterworth'}
filts.bpf = 'eliptic';      % {'none', 'kaiser', 'butterworth', 'eliptic'}
nnl_type  = 'absolute';     % {'absolute', 'rectify'}
% filts.lpf = 'eliptic';      % {'none', 'kaiser', 'butterworth', 'eliptic'}

% downsampling factors
[N, D] = rat(Fs/Fs_new);


%% High-pass filter
switch lower(filts.hpf)
    case 'none'
        % Don't apply high-pass filter
        x_hpf = x;
        
    case 'butterworth'
        error('Not implemented yet!');
        
    case 'kaiser'
        % Butterworth High-pass filter
        fc = 300;           % (Hz) cut off frequency
        fn = Fs/2;          % (Hz) Nyquivst frequency = sample frequency/2;
        order = 4;          % n'th order filter, high pass
        [b_hpf, a_hpf] = butter(order, (fc/fn), 'high');
        x_hpf = filtfilt(b_hpf, a_hpf, x);
        if ~isempty(fignum)
            figure(fignum);
            fignum = fignum + 5;
            fvtool(b_hpf, a_hpf);
        end
    
    otherwise
        error('--> Unrecognized <%s> filter!', filts.hpf);
end
    
    
%% Band-pass filter
switch lower(filts.bpf)
    %{
    case 'kaiser'
        % Belitski et al., 2008
        
        % Kaiser Band-pass filter
        fcuts = [275 400 3000 3500];
        mags  = [0 1 0];
        devs  = [0.001 0.01 0.001];

        [n, Wn, beta, ftype] = kaiserord(fcuts,mags,devs, Fs);
        n = n + rem(n,2);
        b_bpf = fir1(n, Wn, ftype, kaiser(n+1,beta), 'noscale');
        a_bpf = 1;        
        if ~isempty(fignum)
            figure;
            fvtool(b_bpf, a_bpf);
        end

    case 'butterworth'
        fcutlow=100;        %low cut frequency in Hz
        fcuthigh=1700;      %high cut frequency in Hz
        [b_bpf, a_bpf] = butter(order, [fcutlow,fcuthigh]/(fs/2), 'bandpass');
        if ~isempty(fignum)
            figure;
            fvtool(b_bpf, a_bpf);
        end
    %}
    case 'eliptic'
        % filtfilt eliminates the nonlinear phase distortion of an IIR 
        % filter (see: https://www.mathworks.com/help/signal/ug/iir-filter-design.html)
        d_hpf = designfilt('bandpassiir',...
            'FilterOrder', 24, ...
            'PassbandFrequency1', 0.325,...     kHz
            'PassbandFrequency2', 3.0, ...      kHz
            'PassbandRipple', 0.01, ...         dB
            'StopbandAttenuation1', 60,...      dB
            'StopbandAttenuation2', 60, ...     dB
            'SampleRate', 100 );              % kHz 
        
        x_bpf = d_hpf.filtfilt(x_hpf);
        
    otherwise
        error('--> Unrecognized <%s> filter!', filts.bpf);
end


if ~isempty(fignum)
    figure(fignum);
    fignum = fignum + 5;
    clf;
    f = 1e-3*linspace(0, Fs, length(x_bpf))';
    % plot([x_hpf, x_bpf]);
    plot(f, db(abs(fft([x_hpf, x_bpf]))));
    legend('$x_{hpf}$', '$x_{bpf}$');
    xlabel('Frequency (kHz)');
    ylabel('dB');
end



%% Nonlinear operation
% Absolute value; Schnupp et al., 2015.
switch lower(nnl_type)
    case 'absolute'     % Full-wave rectification
        x_nnl = abs(x_bpf);
        
    case 'rectify'      % Half-wave rectification
        x_nnl = max(0, x_bpf);
   
    otherwise
        error('--> Unrecognized <%s> operation!', nnl_type);
end
   
if ~isempty(fignum)
    figure(fignum);
    fignum = fignum + 5;
    clf;
    % plot([x, x_bpf, x_nnl]);
    plot(f, db(abs(fft([x, x_bpf, x_nnl]))));
    legend('$x$', '$x_{bpf}$', '$x_{nnl}$');
    xlabel('Frequency (kHz)');
    ylabel('dB');
end



%% LPF 
d_lpf = designfilt('lowpassiir', 'FilterOrder', 8, ...
         'PassbandFrequency', 0.250,...     % (kHz)
         'PassbandRipple',0.01, ...         % (dB)
         'StopbandAttenuation', 60, ...     % (dB)
         'SampleRate', 100);                % (kHz)

x_lpf = d_lpf.filtfilt(x_nnl);

if fignum
    figure(fignum);
    fignum = fignum + 5;
    clf;
    f = 1e-3*linspace(0, Fs, length(x_bpf))';
    % plot([x_nnl, x_lpf]);
    plot(f, db(abs(fft([x_nnl, x_lpf]))));
    legend('$x_{nnl}$', '$x_{lpf}$');
    xlabel('Frequency (kHz)');
    ylabel('dB');
end




%% Downsample
Fs_dwn = Fs*(D/N);
assert(norm(Fs_dwn-Fs_new)<=5*eps);

% Downsample from Fs to Fs_new
y = resample(x_lpf, D, N);

dt_dwn = 1/Fs_new;      % (sec)
t_dwn = linspace(0, dt_dwn*length(y), length(y));    


if fignum    
    dt = 1/Fs;              % (sec)
    t = linspace(0, dt*length(x_lpf), length(x_lpf));

    figure(fignum);
    clf;    
    plot(t, x_lpf, '-', t_dwn, y, '-');
    xlabel('Time (sec)');
    ylabel('Amp');
    legend('$x_{lpf}$', '$y$');
    
    
    %{
    % f-axis for Fs
    f = 1e-3*linspace(0, Fs, length(x_bpf))';
    
    % f-axis downsampled
    f_dwn = linspace(0, f(end)/N, length(x_bpf)/N)'; 
    
    figure(fignum);
    fignum = fignum + 5;
    clf;    
    subplot(1,2,1);
    plot(f, db(abs(fft(x_lpf))));
    xlabel('Frequency (kHz)');
    ylabel('dB');    
    subplot(1,2,2);
    plot(t, x_lpf);
    legend('$x_{lpf}$');
    xlabel('Time (sec)');
    ylabel('Amp');
    
    figure(fignum);
    fignum = fignum + 5;    
    clf;    
    subplot(1,2,1);
    plot(f_dwn, db(abs(fft(y))));
    xlabel('Frequency (kHz)');
    ylabel('dB');
    subplot(1,2,2);
    plot(t_dwn, y)
    legend('$y$');
    xlabel('Time (sec)');
    ylabel('Amp');
    %}
    

end

























