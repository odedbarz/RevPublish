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
addRequired(p, 'type', @(x) ischar(x) || isstring(x));            

% Optional:
% - CC:
addParameter(p, 'fn', '', @(x) ischar(x) || isstring(x));

% - SVD
addParameter(p, 'Y', [], @isnumeric);            
addParameter(p, 'n_svd', 1, @isnumeric);         % # of singular values (SVs) to use

% - RND & SPK & NOSPK
addParameter(p, 'N', [], @isnumeric);           % total number of units    

% - SPK & NOSPK
% addParameter(p, 'tbl_data', [], @istable);           % total number of units    
addParameter(p, 'fn_template', [], @iscell);           % total number of units    


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
        [~, sorted_list] = sort(tbl_BF.BF_cc, 'descend');
        
        varargout{1} = tbl_BF;
        
        
    case 'FILE'
        fpath = load.path_to_data('Analysis');
        [~, fname, ~] = fileparts(pars.fn); 
        fn    = fullfile(fpath, [fname, '.mat']);
        assert(~isempty(dir( fn )), '--> ERROR: can''t find the DATA file!');
        dummy = load( fn );
        sorted_list = dummy.unit_list;

        
        
    case 'SVD'
        % Orthogonalize the measurements
        %[U, S, ~] = svds(zca(pars.Y), pars.n_svd);
        [U, S, ~] = svds(pars.Y, pars.n_svd);
        Z = U(:, 1:pars.n_svd); % * S(1:pars.n_svd, 1:pars.n_svd);

        % Re-arrange the units; best units to contribute the reconstruction at the
        % beginning
        P = pars.Y'*Z;      % projections
        [~, sorted_list] = sort(P, 'descend');

        % flatten the 2D matrix
        sorted_list = unique(sorted_list', 'stable')';

        % Singular values
        varargout{1} = diag(S(1:pars.n_svd, 1:pars.n_svd));
        
        
    case 'RND' 
        if isfield(pars, 'N') && ~isempty(pars.N)
            N = pars.N;             
        elseif isfield(pars, 'Y')
            N = size(pars.Y,2); 
        else
            error('Please enter number of units (N)!');
        end
        
        sorted_list = randperm(N);
        varargout{1} = {};
        
        
    case 'SPK'     
        data_type = pars.fn_template{3};
        fn_path = pars.fn_template{1};
        
        fn_su = sprintf(pars.fn_template{2}, 'SU');
        fn_SU_full = fullfile(fn_path, fn_su);
        data = load(fn_SU_full);
        tbl_su = data.(sprintf('tbl_%s', 'SU'));
        
        fn_mua = sprintf(pars.fn_template{2}, 'MUA');
        fn_MUA_full = fullfile(fn_path, fn_mua);
        data = load(fn_MUA_full);
        tbl_mua = data.(sprintf('tbl_%s', 'MUA'));
        
        plausible_neurons = intersect(tbl_su.neuron, tbl_mua.neuron);
        
        % Get the correct lines in the relevant table
        tbl_slc = eval(sprintf('tbl_%s', lower(data_type)));   % tbl_mua OR tbl_su
        plausible_units = arrayfun(@(X) find(tbl_slc.neuron == X), plausible_neurons);      
        
%         if isfield(pars, 'N') && ~isempty(pars.N)
%             N = pars.N;             
%         elseif isfield(pars, 'Y')
%             N = length(plausible_units); 
%         else
%             error('Please enter number of units (N)!');
%         end
        
        % Get N random units (if require)
        %ind_units = randperm( length(plausible_units), N );
        %sorted_list = plausible_units(ind_units);
        sorted_list = plausible_units;
        arargout{1} = {};

        
        
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




