function S = Impale_struct(filename)
%
%   function S = load.Impale_struct(filename)
%
% Input:
%   filename: (str) filename to Impale's data structure.
%
% Output:
%   S: Impale's structure
%
% Description:
% Sets the tabs for all selected measurements.
%

import load.*

% Load Impale structures
warning off
S = load( filename );
warning on

try 
    warning('Off');
    S = load( filename );
    warning('On');    
catch
    o.err('--> ERROR in [load_spikes_from_Impale.m]:\n\t Could not load this file <%s> !!\n', filename);
end


assert(isstruct(S) && isfield(S, 't') && isfield(S, 'ch'),...
    '--> ERROR in [load_spikes_from_Impale.m]:\n\t <filename> MUST point into Impale''s structure!!\n')






