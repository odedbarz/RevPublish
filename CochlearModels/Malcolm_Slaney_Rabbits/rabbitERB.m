function [rabbitERBres, Qerb] = rabbitERB(f)
%Rabbit ERBs 
%
%This function finds rabbit ERBs based on the calculations of Carney,
%Zhang, and Glasberg and Moore (BASED ON GAI et al. 2007) but with a new
%way of fitting the Borg, 1988 et al. data 
orig_f=f;

f=f/1000; %convert to kHz

%Q10 from fit of Borg et al. Q10 data for Rabbit AN tuning curves
Q10=3.62*(f.^3.65./(f.^3.65+3.27^3.65))+2.022;

%Based on empirical relationship between Q10 and Q_ERB (footnote 6 in
%Shera & Guinan 2003)
Qerb=1.75*Q10;

%Calculate ERB bandwidth.
rabbitERBres = orig_f./Qerb;


