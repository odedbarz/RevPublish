function [xSD] = sdgn(x, n_win)
% 
%
% Description:
% Synaptic depression Gain normalization
%
% Mesgarani et al., 2014, Mechanisms of noise robust representation of
% speech.

%% Synaptic depression (SD) model
%{
% x = rtr;
% y = rest;
% a = lpc(x, 10);
% % est_y = filter([0 -a(2:end)], 1, y);
% est_y = filter([0 -a(2:end)], 1, y);
% est_y = circshift(est_y, -1);

% x = rts;
% y = ones(size(rest));   y(701:900) = y(701:900) + hann(200); 'DEBUG'
y = rest;
n_win = 150;
W = hann(n_win)'; 
W = W((end+1)/2:end);
W = W/sum(W);
% W = W/norm(W)^2;
a = W(1:end-1);
Vth = 1;
At    = Vth * filtfilt([0, a(2:end)], 1, y);
% At    = Vth * filter([0, 1/2], 1, y);
est_y = y - At;
% est_y = filtfilt(a(2:end), 1, y);


fprintf('\n-> ARMA model\n');
cc = corrcoef(rts0, rest);
fprintf('--> corrcoef(rts0, rest): %.3f\n', cc(1,2));
cc = corrcoef(rts0, est_y);
fprintf('--> corrcoef(rts0, est_y): %.3f\n', cc(1,2));

figure(14);
clf;
plot([y, At, est_y]);
% plot([rts0,  rest, est_y]);
% plot(zca([rts0, rts, rest, est_y]));
%}


%%
sd_type = 3
switch sd_type
    case 1  % lpc
        [Az, G] = lpc(x, 13);
        Bz = -Az(:, 2:end);     % B(z) = 1 - A(z)
        
    case 2  % window
        n_win = 2*fix(n_win/2) + 1; % always an odd number
        W = hanning(n_win)';        % hanning() doesn't have zeros at its margins like hann()
        W = W((end+1)/2:end);
        W = W/sum(W);
        Bz = W;                    % drop the zero at the end
    case 3
        n_win = 2*fix(n_win/2) + 1; % always an odd number
        W = hanning(n_win)';        % hanning() doesn't have zeros at its margins like hann()
        W = W((end+1)/2:end);
        W = W/sum(W);
        Bz = [0, W(2:end)];
        
end

    
xSD = nan(size(x));
for k = 1:size(x,2)
    xSD(:,k) = filter(Bz(k,:), 1, x(:,k));
end

% x   = x(2:end,:);
% xSD = xSD(1:end-1,:);
% nt  = size(x,1);


















