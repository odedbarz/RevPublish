function ynew = PwLPF(y)

ref_amp = 0.01*min(min(y(end/2:end,:)));
N = 2*100;   % max window size

W = hann(N);
W = W(end/2:end);
W = W/sum(W);
W = [W; zeros(size(W))];

ylpf = fftfilt(W, y);

ynew = max((y-ylpf)./(y + 1e-6), ref_amp);

