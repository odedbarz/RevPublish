function A = doubly_centered_dist(x)
%
% Calculate doubly centered distance matrix using MATLAB's vectorization 
%
    a = pdist2(x,x);
    ajbar = mean(a);
    akbar = mean(a,2);
    abar = mean(mean(a))*ones(size(a));
    % ajbar = ones(size(mrow))*mcol;
    % akbar = mrow*ones(size(mcol));
    % abar = mean(mean(a))*ones(size(a));
    % A = a - ajbar - akbar + abar;
    A = a - ajbar - akbar + abar;

end