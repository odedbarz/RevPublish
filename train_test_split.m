function [X_train, X_test, y_train, y_test, split_st] = train_test_split(X, y, varargin)
%
%   function [X_train, X_test, y_train, y_test, splits] = train_test_split(X, y, varargin)
%
% [n_splits]    : (1x1 int) number of splits to split the data. NOte that 
%                 this value is ignored if SPLIT_TIMES is given.
%
% [split_time_idx] : (T x 2) an array of starting (first column) and ending
%                 (second column) of the chunks (time intervals) to split.
%                 If SPLIT_TIMES is given (not empty), then N_SPLITS is override
%                 with the number of rows in SPLIT_TIMES.
%   
% [test_grp]    : (1x1 int) the number of test group to use out of all
%                 N_SPLITS groups
% [dim_x]       : (1x1 int; default=first non-singleton dimension of X) splits 
%                 the matrix x along this dimension. 
% [dim_y]       : (1x1 int; default=first non-singleton dimension of y) splits 
%                 the matrix y along this dimension. 
%
% This implementation follows the same function from Python, see
%   sklearn.model_selection.train_test_split.
%
%


%% Check the Input
p = inputParser;

addRequired(p, 'X', @ismatrix);
addRequired(p, 'y', @ismatrix);

% # frequencies
addOptional(p, 'test_grp', 12, @(x) isscalar(x));    % 
addOptional(p, 'n_splits', 12, @(x) isscalar(x));
addOptional(p, 'split_time_idx', [], @(x) isnumeric(x));
addOptional(p, 'dim_x', 2, @(x) isnumeric(x));
addOptional(p, 'dim_y', 1, @(x) isnumeric(x));

% 
addOptional(p, 'verbose', 0, @(x) isnumeric(x) || islogic(x));

parse(p, X, y, varargin{:});
pars = p.Results;


% Check the data
len_x = size(X, pars.dim_x);
len_y = size(y, pars.dim_y);
assert(len_x == len_y, '--> ERROR: the provided data have different length and cannot be split along these dimensions!');


%% Split the data
if isempty(pars.split_time_idx)
    % Group indices; split into (approximately) same size N_SPLITS groups;
    % Deviation could be of about N_SPLITS samples at the last group
    split_st.idx = fix((0:(len_x-1))/ceil(len_x/pars.n_splits)) + 1;

else
    % In this case, each interval can have different length. So, set each
    % interval index (k) separately.
    pars.n_splits = size(pars.split_time_idx, 1);
    split_st.idx = zeros(1, len_x);
    for k = 1:pars.n_splits
        chunk_k = pars.split_time_idx(k,1):pars.split_time_idx(k,2);        
        split_st.idx(chunk_k) = k;
    end    
    
end

split_st.n_grps  = pars.n_splits;
split_st.test    = pars.test_grp;
split_st.n_splits= pars.n_splits;
split_st.train   = setdiff(1:split_st.n_splits, split_st.test);


if pars.verbose
    a = arrayfun(@(x) nnz(split_st.idx==x), 1:split_st.n_splits);
    fprintf('--> n_splits           : %d (# of groups)\n', split_st.n_splits);
    fprintf('--> # samples by group :\n\t--> [%s]\n', num2str(a, ' %d'));
    fprintf('--> max # of deviations: %d samples\n', max(median(a) - min(a)));
    fprintf('\n');
end


%% Check the input: remove empty columns in y 
% (this measurement might not include the requested DRR case)
is_empty_col = (0 == sum(y)) | ( isnan(sum(y)) );
y = y(:,~is_empty_col);
if pars.verbose
    if 0 < nnz(is_empty_col)
        warning('[train_test_split.m]: One or more measurements might not include the requested DRR case!');
    end
end


%% Split the data
% Always perform the split along the first dimension
if 1 < pars.dim_x 
    X = X';
end

if 1 < pars.dim_y 
    y = y';
end

X_train = X(split_st.idx ~= pars.test_grp, :);
X_test  = X(split_st.idx == pars.test_grp, :);
y_train = y(split_st.idx ~= pars.test_grp, :);
y_test  = y(split_st.idx == pars.test_grp, :);

% switch back, if necessary
if 1 < pars.dim_x 
    X_train = X_train';
    X_test  = X_test';
end

if 1 < pars.dim_y 
    y_train = y_train';
    y_test  = y_test';
end




