function [Yenv, Yfilt] = LFP_filtbank(Y, fs, varargin)
%
%   function [Yenv, Yfilt] = LFP_filtbank(Y, fs, [fc], [bw], [N])
%
% Description:
%


%% Set the inputs
p = inputParser;

addRequired(p, 'Y', @isnumeric);            
addRequired(p, 'fs', @isnumeric);      % (Hz) samplnig rate of the input stimulus <y>

% Default frequencies (Stark & Abeles, 2007)
%     alpha       : 1-13 Hz     -> ~6 Hz, bw=1/3
%     beta        : 13-30 Hz    -> ~16 Hz, bw=1/3
%     gamma       : 30-60 Hz    -> ~45 Hz, bw=1
%     higher gamma: 60-100 Hz   -> ~75 Hz, bw=1  
addOptional(p, 'fc', [  6,  16, 45, 75], @isnumeric);       
addOptional(p, 'bw', [1/3, 1/3,  1,  1], @isnumeric);       
addOptional(p, 'N', 3, @isnumeric);             
addOptional(p, 'fignum', [], @isnumeric);             

parse(p, Y, fs, varargin{:});

pars  = p.Results;
fignum= p.Results.fignum;    


%%
[n_smp, n_y] = size(Y);
n_filt= length(pars.fc);
Yenv  = nan(n_smp, n_filt, n_y);
Yfilt = nan(n_smp, n_filt, n_y);

for k = 1:n_y
    [Yenv(:,:,k), Yfilt(:,:,k)] = octave.filter_bank(Y(:,k), fs, pars.fc, pars.bw, pars.N);
end


%% Plot
if isempty(fignum)
    return;
end

figure(fignum);
clf;
[nx, ny, ~] = size(Yfilt);
offset = 1.5*max(max(abs(zca(Yfilt(:,:,1)))));
yoff = offset * ones(nx,ny) * diag([0:ny-1]);
% plot(zca(Yfilt(:,:,1)) + yoff);
% hold on
plot(zca(Yenv(:,:,1)) + yoff);
% hold off
set(gca, 'YTick', offset * [0:ny-1]);
set(gca, 'YTickLabel', {'$\alpha$', '$\beta$', '$\gamma$', '$\gamma_{high}$'});
















