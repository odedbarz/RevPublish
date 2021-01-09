function [sorted_list, varargout] = find_best_unit_set(type, varargin)
%
%   function sorted_list = find_best_unit_set(Y, varargin)
%
% Description:
%   Orthogonalize the measurements to remove multicolinearity. Then project
%   the orthogonalized measurements (P) on the orthogonal axes to get a
%   set of units that best contributes to the spectrogram's reconstruction.
%

%% Set the inputs
p = inputParser;

% Required
addRequired(p, 'type', @(x) isstr(x) || isstring(x));            

% Optional:
% - CC:
addParameter(p, 'fn', '', @(x) isstr(x) || isstring(x));

% - SVD
addParameter(p, 'Y', [], @isnumeric);            
addParameter(p, 'n_svd', 1, @isnumeric);         % # of singular values (SVs) to use

% - Verbose
addParameter(p, 'fignum', [], @isnumeric);           
 
parse(p, type, varargin{:});
pars = p.Results;


 

%%
switch upper(pars.type)
    case 'CC' 
        [fpath, fname, fext] = fileparts(pars.fn);        
        assert(~isempty(dir( pars.fn )), '--> ERROR: can''t find the DATA file!');
        fn_bf = fullfile(fpath, [fname, '_BF', fext]);
        assert(~isempty(dir( fn_bf )),...
            '--> ERROR: can''t find the _BF file! If needed, create the file with best_envelope_frequency.m');
        dummy = load( fn_bf );
        tbl_BF = dummy.tbl_BF;
        [~, sorted_list] = sort(tbl_BF.BF_cc);
        
        varargout{1} = tbl_BF;
        
        
    case 'SVD'
        % Orthogonalize the measurements
        [U, S, ~] = svds(zca(Y), pars.n_svd);
        Z = U(:, 1:pars.n_svd); % * S(1:pars.n_svd, 1:pars.n_svd);

        % Re-arrange the units; best units to contribute the reconstruction at the
        % beginning
        P = Y'*Z;      % projections
        [~, sorted_list] = sort(P, 'descend');

        % flatten the 2D matrix
        sorted_list = unique(sorted_list', 'stable')';

        % Singular values
        sv = diag(S(1:pars.n_svd, 1:pars.n_svd));
        

    otherwise
        error('ERROR: type <%s> is unrecognized!', pars.type);
        
end



%% Plot\Debug
if ~isempty(pars.fignum)    
    figure(pars.fignum);
    clf;
    plot(P, '.-');
    xlabel('Unit Number');
    ylabel('$Y^T\cdot Z$');
    title('Projections over Orthogonal Axes');
end




