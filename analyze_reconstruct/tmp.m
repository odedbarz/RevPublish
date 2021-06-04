function ynew = win(y)

ref_amp = 0.01*min(y(end/2:end));
N = 2*50;   % max window size

W = hann(N);
W = W(end/2:end);

ylpf = conv(y, W);

y = max((y-ylpf)/(y + 1e-6), ref_amp);

