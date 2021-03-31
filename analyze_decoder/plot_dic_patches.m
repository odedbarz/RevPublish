function [ax, patches_on_x, patches_on_y] = plot_dic_patches(D, m, n, dt, faxis)
%
% function [ax, patches_on_x, patches_on_y] = plot_dic_patches(D, m, n)
%
% Description;
%   Displays patches (columns of D). The patches are sorted along the
%   columns, in the usual column-first "Matlab's way".
%
% Oded Barzelay, 12/18/2019

%% Set the inputs
[L, dic_dim2] = size(D);

% Limit the number of features to show
MAX_FEATURES = 64;
n_patch = min(dic_dim2, MAX_FEATURES);
patches_on_x = ceil(sqrt(n_patch))-1;
patches_on_y = patches_on_x;
n_patch      = patches_on_x * patches_on_y;
patch_num   = randperm(dic_dim2, n_patch);
patch_num   = sort(patch_num);
% D = D(:, patch_num);

if 3 > nargin || isempty(n)
    n = L/m;
    assert(n == fix(L/m), '--ERROR: L must be a product of two INTEGERS, m & n!');
end

assert( L == m*n, '--> ERROR:  L != m*n; The dimension of the patches disagree with the number of rows in D! ');

if 4 > nargin
    dt = [];
end

if 5 > nargin
    faxis = [];
end


%%
% # of patches along the x & y dimensions

% Prepare the figure 
gcf;
clf;

% Holder for the axes
ax = zeros(patches_on_x, patches_on_y);


for k = 1:(patches_on_x * patches_on_y)
    patch_k = reshape(D(:,patch_num(k)), m, n);
    
    % Normalize the patch before plotting it
    patch_k = patch_k/(eps + max(abs(patch_k(:))));
        
    ax(k) = subplot(patches_on_x, patches_on_y, k);
    imagesc(patch_k);
    set(gca, 'YTickLabel', '');
    set(gca, 'XTickLabel', '');
    set(gca, 'YDir', 'normal');
    set(gca, 'Box', 'off');
end

%ax = ax';   % to match between the display and the layout of the matrix

fontsize = 11;

% Set scales for the 
if ~isempty(faxis)
    axes( ax(1,end) );
    set(gca, 'YTick', [1, m]);
    yticklabel = arrayfun(@(D) sprintf('%g', 1e-3*D), faxis([1, end]), 'UniformOutput', 0);
    set(gca, 'YTickLabel', yticklabel)
    ylabel('Frequency (kHz)');
    set(gca, 'FontSize', fontsize);
end

if ~isempty(dt)
    set(gca, 'XTick', [1, n])
    xticklabel = arrayfun(@(D) sprintf('%g', D), dt*([1, n]), 'UniformOutput', 0);
    set(gca, 'XTickLabel', xticklabel)
    xlabel('Lags (ms)');
    set(gca, 'FontSize', fontsize);
end



