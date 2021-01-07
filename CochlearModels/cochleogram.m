%
% cochleogram.m
%
function cochleogram(yt, Fs, t0)
%
%   function cochleogram(yt, Fs, t0, t1)
%

% t0 = 0;

%% 
% summation_window              : 0.025 sec
sumwin = 0.010;

% hop between successive windows: 0.010 sec 
hopdist = 0.005;

% number of gammatone channels  : 64
Ngamma = 64;

% lowest frequency              : 50 Hz
lowfrew = 50;   

% highest frequency             : fS/2
highfreq = Fs/2;

[D2, cfArray] = gammatonegram(yt, Fs, sumwin, hopdist, Ngamma, lowfrew, highfreq, 0);
% imagesc(20*log10(D2));


%% Plot

clf;
dummy = (1:8)+10*[0:3]';
subplot(5,10,dummy(:));

% D2plot = 20*log10(abs(D2));
D2plot = D2;

tt = t0 + (0:1/Fs:(size(yt,1)-1)/Fs);  % time axis
% ff = linspace(0, Fs-1/length(tt), length(tt));  % frequency domain
ttt = linspace(t0, max(tt)-1/Fs, size(D2,2));
fff = log10( cfArray );

[TTT,FFF] = meshgrid(ttt, fff);
surface(TTT,FFF,D2plot,'EdgeColor','none');     
ax(1) = gca;
ax(1).XTickLabel = '';
ax(1).YScale = 'linear';
ax(1).YDir = 'normal';      % 'reverse'
axis tight
ax(1).YTickLabel = num2str(1e-3*10.^(str2double(ax(1).YTickLabel)),'%.2f');
set(gca, 'fontsize', 16);
% xlabel('Time (s)');
ylabel('CF (kHz)', 'FontSize', 20);
ylim1 = ylim;
aux.ctitle('Cochleogram $C(t,f)$' , ...
    sprintf('(Window: %g ms, Overlap %g ms, Gammatones: %d)', 1e3*sumwin, 1e3*(sumwin-hopdist), Ngamma) ...
    );
% caxis([0, 0.3])

dummy = (9:10)+10*[0:3]';
subplot(5,10,dummy(:));
plot(sum(D2plot,2), fff);
ax(2) = gca;
ax(2).YTickLabel = '';
%ax(2).YTickLabel = num2str(1e-3*10.^(str2double(ax(2).YTickLabel)),'%.2f');
axis tight
ax(2).YScale = 'linear';
ax(2).YDir = 'normal';
ylim(ylim1);
linkaxes(ax, 'y')
ylabel('');
title('$C(f)$')

subplot(5,10,(1:8) + 10*4);
plot(ttt, sum(D2plot,1));
ax(2) = gca;
% ax(2).YTickLabel = '';
%ax(2).YTickLabel = num2str(1e-3*10.^(str2double(ax(2).YTickLabel)),'%.2f');
axis tight
ax(2).YScale = 'linear';
ax(2).YDir = 'normal';
% linkaxes(ax, 'x')
xlabel('Time (s)', 'FontSize', 20);
ylabel('$C(t)$')










