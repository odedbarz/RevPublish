%
% test_sdgn
%

rng(1)

% nt = 2000;
% xr = randn(nt,1);
% b = fir2(128, [0 0.05 0.1 1], [1 1 0 0]);
% x = filtfilt(b, 1, xr);

x = H(:,:,2);


%% xSD
% xsd = sdgn(x, n_win);
x0 = x(:,3);
x_ = x(:,1);

n_win = 550;

% n_win = 2*fix(n_win/2) + 1; % always an odd number
% W = hanning(n_win)';        % hanning() doesn't have zeros at its margins like hann()
% W = W((end+1)/2:end);
% W = W/sum(W);
% Bz = -W(2:end);                    % drop the zero at the end
% xsd = x_ + filter(Bz, 1, x_);

n_win = 2*fix(n_win/2) + 1; % always an odd number
W = hanning(n_win)';        % hanning() doesn't have zeros at its margins like hann()
W = W((end+1)/2:end);
W = W/sum(W);
Bz = [0, W(2:end)];
xsd = x_ - filter(Bz, 1, x_);
xsd = max(0, xsd);

%
figure(11);
clf
plot([x0, x_, xsd]);

corrcoef(x0, xsd)
corrcoef(x_, xsd)






