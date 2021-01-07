function [Nneu, Nses] = how_many_neurons(T)
%
%   function [Nneu, Nses] = how_much_neurons(T)
%
% Input:
%   T: A table that list all measurements.
% 
% Output:
%   N: # of neurons (** not measurements **) in T.
%
%

sessions = T{:,'session'};      % # of sessions
units    = T{:,'unit'};         % # of units
spikechan= T{:,'spikechan'};         % # of units
Nneu = length( unique(1e4*sessions + 1e2*units + spikechan) );

Nses = length(unique(sessions));



