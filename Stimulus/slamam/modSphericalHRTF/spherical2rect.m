function [x,y,z] = spherical2rect(r,az,el)

%[x,y,z] = spherical2rect(r,az,el)
%
%Converts spherical coordinates (r,azimuth,elevation) to rectangular
%coordiantes in (x,y,z) space
%
%

%step 1: convert angles to radians
az = pi*az/180;
el = pi*el/180;

az
el

x = r * cos(az) * sin(el);
y = r * sin(az) * sin(el);
z = r * cos(el);
