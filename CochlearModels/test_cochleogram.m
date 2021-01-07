%
% cochleogram.m
%

clc

%% Load the data
fn = 'C74 cries 2.wav';
% fn = 'call_1.mp3';

fnfull = ['.data\Rabbit Calls/', fn];

[yt, Fs] = audioread(fnfull);   
yt = yt(:,1);


%%
figure(11);

t0 = 6;   % s
t1 = 7; 

yt1 = yt([t0*Fs:t1*Fs]);
cochleogram(yt1, Fs, t0);

% t0 = 6615/Fs;   % s
% yt1 = yt([6615:19845]);
% cochleogram(yt1, Fs, t0);



%% fft
% %{
tt = t0 + (0:1/Fs:(size(yt,1)-1)/Fs)';  % time axis
ff = linspace(0, Fs-1/length(tt), length(tt));  % frequency domain

Yf = fft(yt);

figure(51);
clf;
subplot(3,1,1);
plot(tt, yt);
% hold on
% plot(tt, abs(hilbert(yt)), '-');
% hold off
xlabel('Time (sec)');
ylabel('y(t)');
axis tight
title(mName2latex( sprintf('Fs: %g Hz (file: %s)', Fs, fn) ));

subplot(3,1,2);
plot(1e-3*ff, abs(Yf));
xlabel('Frequency (kHz)');
ylabel('$Y(f)$');
xlim([0, 1e-3*Fs/2]);

subplot(3,1,3);
plot(1e-3*ff, 20*log10(abs(Yf)));
xlabel('Frequency (kHz)');
ylabel('$20\cdot log_{10}\{Y(f)\}$');
xlim([0, 1e-3*Fs/2]);
%}



%%
figure(15);
clf;

Ed = emd(yt1)';

imagesc( gammatonegram(yt1, Fs, 0.005, 0.001, 2*64, 50, Fs/2, 0) )
hold on
contour( gammatonegram(Ed(:,1), Fs, 0.005, 0.001, 2*64, 50, Fs/2, 0), 3, 'Color', rpalette('new02') )
contour( gammatonegram(Ed(:,2), Fs, 0.005, 0.001, 2*64, 50, Fs/2, 0), 3, 'Color', rpalette('new03') )
contour( gammatonegram(Ed(:,4), Fs, 0.005, 0.001, 2*64, 50, Fs/2, 0), 3, 'Color', rpalette('new05') )
hold off













