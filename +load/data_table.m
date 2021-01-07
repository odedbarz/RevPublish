function T = data_table(data_in_st, varargin)
%
%   function T = load.data_table(data_st)
%

import load.*


%% Set the input structure   
% Default input
data_st.path        = path_to_data; 
data_st.fn          = 'LabNoteBook.xlsx';
data_st.operator    = 'OBZ';
data_st.ID          = 'C74';
data_st.measType    = '';
data_st.session     = [];     
data_st.unit        = [];      
data_st.measNum     = [];     
%data_st.syncchan   = [];
%data_st.spikechan  = [];     

% Assign input(s), if given
if exist('data_in_st', 'var') && ~isempty(data_in_st) && isstruct(data_in_st)
    
    % fields_st = fieldnames(data_st);
    fields_in = fieldnames(data_in_st);
    for k = 1:length(fields_in)
        assert(isfield(data_st, fields_in{k}), '--> ERROR in [load.data_table]: The field <%s> is invalid!!!', fields_in{k});
        data_st.(fields_in{k}) = data_in_st.(fields_in{k});
    end
    
end



%%
filename_full = [data_st.path, data_st.fn];

if isempty(dir(filename_full))
    fprintf('--> ERROR in [load.data_table.m]: can''t find the requested file!!!\n');
    fprintf('----> filename: %s\n', filename_full);
    error('can''t find the requested file!!!');
end

[~, ~, data_st.filetype] = fileparts( data_st.fn );

T = readtable(filename_full, 'Sheet', data_st.ID);


%% Sift the data
if ~isempty(data_st.ID)
    idx = cellfun(@(ID) strcmpi(ID, data_st.ID), T.ID, 'UniformOutput', 1);
    T = T(idx,:); 
end

if ~isempty(data_st.measType)
    idx = cellfun(@(TYPE) strcmpi(TYPE, data_st.measType), T.measType, 'UniformOutput', 1);    
    T = T(idx,:); 
end

if ~isempty(data_st.session)
    idx = T.session == data_st.session;
    T = T(idx,:); 
end

if ~isempty(data_st.unit)
    idx = T.unit == data_st.unit;
    T = T(idx,:); 
end

if ~isempty(data_st.measNum)
    idx = T.measNum == data_st.measNum;
    T = T(idx,:); 
end


%% Add ROWs
row = (1:size(T,1))';
T = [table(row), T];



















