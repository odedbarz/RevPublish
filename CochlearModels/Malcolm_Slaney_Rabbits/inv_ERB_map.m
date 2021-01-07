function f = inv_ERB_map(n)
% INV_ERB_MAP - Maps ERB number to frequency in Hz
%  Usage: f = inv_ERB_map(n)
%
f = (10.^(n/21.4) -1)/4.37e-3;