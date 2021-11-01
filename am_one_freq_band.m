function [x_am, abs_xf, f] = am_one_freq_band(x, dt, nf, plot_output)
    
nx = length(x);
if 3 > nargin 
    nf = nx;
end
if 4 > nargin
    plot_output = false;
end

x = x(:);

fs = 1/dt;
f = [0:fs/(nf-1):fs]';
abs_xf = abs(fft(x, nf));    
x_am = abs_xf/abs_xf(1);

if plot_output    
    subplot(1,2,1);
    plot(f, abs_xf); 
    xlim([0, fs/2]);
    title('|FFT(x)|');
    subplot(1,2,2);
    plot(f, x_am); 
    xlim([0, fs/2]);
    xlabel('Frequency (Hz)');
end

