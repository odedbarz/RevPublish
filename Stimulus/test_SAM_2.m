


clc

Fs  = 100e3;                    % (Hz)
Nt  = ceil(Fs*0.5);             % (smp)
t   = linspace(0, 0.5, Nt)';    % (sec)
f   = linspace(0, Fs, Nt)';     % (Hz)
CF  = 4e3;

% Envelope (sine)
fenv = 32;                      % (Hz)
s1 = sin(2*pi*t*32);
S1 = 1/Nt*fft(s1);

% The neural 'spectral receptive field' filter
Hd = BW_BPF(CF);

% SAM
x1 = s1 .* randn(size(s1));
X1 = 1/Nt*fft(x1);

% The SAM as seen by the IC neuron
xg1 = Hd.filter(x1);

figure(11);
clf;
subplot(1,2,1);
plot(1e3*t, [x1, xg1]);
xlim([0, 100]);
xlabel('Time (msec)');
ylabel('Amp.');

subplot(1,2,2);
plot(1e-3*f, 20*log10(abs(fft([x1, xg1]))));
aux.vline(1e-3*CF);
xlim(1e-3*[2e3, 24e3]);
xlabel('Frequency (kHz)');
ylabel('20log_{10}()');
legend('SAM', sprintf('filter(SAM, CF=%g kHz)', 1e-3*CF), 'Location', 'southeast');







