function [optimal_sorted, P, sv] = find_best_unit_set(Y, varargin)
%
%   function [optimal_sorted, P] = find_best_unit_set(H, train_drr, varargin)
%
% Description:
%   Orthogonalize the measurements to remove multicolinearity. Then project
%   the orthogonalized measurements (P) on the orthogonal axes to get a
%   set of units that best contributes to the spectrogram's reconstruction.
%

%% Set the inputs
p = inputParser;

% Required
% addRequired(p, 'obj', @isobject);            
addRequired(p, 'Y', @(x) isnumeric(x));            
%addRequired(p, 'train_drr', @isnumeric);            

% Optional
addOptional(p, 'n_svd', 1, @isnumeric);         % # of singular values (SVs) to use
addOptional(p, 'fignum', [], @isnumeric);           
 
parse(p, Y, varargin{:});
pars = p.Results;


%% Orthogonalize for the best set
% Define training matrix
% Y = squeeze( H(:, train_drr, :) );
% Y = squeeze( H(:, train_drr, :) );

% Orthogonalize the measurements
[U, S, ~] = svds(zca(Y), pars.n_svd);
Z = U(:, 1:pars.n_svd); % * S(1:pars.n_svd, 1:pars.n_svd);

% Re-arrange the units; best units to contribute the reconstruction at the
% beginning
P = Y'*Z;      % projections
[~, optimal_sorted] = sort(P, 'descend');

% flatten the 2D matrix
optimal_sorted = unique(optimal_sorted', 'stable')';

% Singular values
sv = diag(S(1:pars.n_svd, 1:pars.n_svd));

%% Plot\Debug
if ~isempty(pars.fignum)    
    figure(pars.fignum);
    clf;
    plot(P, '.-');
    xlabel('Unit Number');
    ylabel('$Y^T\cdot Z$');
    title('Projections over Orthogonal Axes');
end




