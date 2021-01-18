function ent = spectral_flatness(Sft, dim, eps)
%
% function ent = spectral_flatness(Sft)
%
% Description:
% Flatness entropy, also know as Wiener entropy, is a measure of how "tone-
% like" a spectrum is compared to "noise-like" spectrum.

if nargin < 2 || isempty(dim)
    dim = 1;
end

if nargin < 3 || isempty(eps)
    eps = 1e-12;
end


% num = geomean(Sft, dim);
num = exp(mean(Sft, dim));
den = mean(exp(Sft), dim);  % Sft is already on LOG scale
ent = num./(den + eps);
