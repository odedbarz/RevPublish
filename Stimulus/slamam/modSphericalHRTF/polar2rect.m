function [x, y, z] =  polar2rect(r, theta)

% elevation is always 0 degrees
z = 0;

y = r * cos(pi*theta/180);
x = r * sin(pi*theta/180);



