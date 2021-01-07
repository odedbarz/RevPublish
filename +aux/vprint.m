function vprint(verbose, varargin)
%
%   function vprint(verbose, varargin)
%
% Description:
% Conditional print: print only when the 'verbose' flag is ON.

if verbose 
    fprintf( varargin{:} ); 
end
