function [valid_file_names, valid_files] = scan_file_template(file_path, template, verbose)
%
%   function [valid_file_names, valid_files] = scan_file_template(file_path, template, verbose)
%

if 3 > nargin
    verbose = 0;
end

% Scan for files to use
file_list       = dir(file_path);
n_files         = size(file_list, 1);
assert(n_files>0, '--< ERROR: can''t find files in the given path!');
valid_files     = false(1, n_files);
parameter_list  = fields(template);
n_pars          = length(parameter_list);
valid_file_names= {};

for k = 1:length(file_list)   
    fn = file_list(k).name;
    
    % Check that the k'th file name contains all parameters    
    % * TEMPLATE must contain a "prefix" field
    valid_file_k = contains(fn, sprintf('%s', template.prefix));
    for m = 2:n_pars
        if ~valid_file_k
            break;
        end        
        check_m = contains(fn, sprintf('%s(%d)', parameter_list{m}, template.(parameter_list{m})));
        valid_file_k = valid_file_k & check_m;
        if verbose
            aux.vprint(~valid_file_k, '--> [scan_file_template.m]: can''t find %s\n', parameter_list{m});
        end
    end
    valid_files(k) = valid_file_k;    
    
    if valid_files(k)
        valid_file_names{end+1} = fn;
    end
end
assert(0<nnz(valid_files), '--> ERROR: couldn''t find ANY file in this path!')

if verbose
    % Total # of files found to use
    n_files = nnz(file_idx);
    fprintf('--> Found %d data files\n', n_files);
end


