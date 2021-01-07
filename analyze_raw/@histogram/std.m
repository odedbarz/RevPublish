function h = std(h)
% HISTOGRAM/STD - Standard deviaton of histogram 
% STD(H) returns a histogram whose data are the standard deviation of those in H
%
h.data = full(std(h.data));
