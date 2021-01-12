clc

nn = 30
len = 600
len_ms = binwidth * len
I = nn+(1:len);

gof = goodness(X_test(:,I), X_est(:,I))

figure(1)
imagesc([X_test(:,I); X_est(:,I)])
title(sprintf('CC: %g', gof.CC))


fidx = 25
fprintf('-> f(n): %g Hz\n', f(fidx));

figure(2)
plot(zca([X_test(fidx,I)', X_est(fidx,I)']));
title(sprintf('CC: %g', gof.CC))




%%

Sdry = spec_st.Sft{3};
Sdrr = spec_st.Sft{4};
ydry = Sdry(fidx,:)';
ydrr = Sdrr(fidx,:)';


len_ydry = length(ydry);
div = 50;
assert( (len_ydry/div) == round(len_ydry/div), '-- Choose a divisable number' )
v = 1:(len_ydry/div):len_ydry;
vidx = v(9) + (1:diff(v(1:2)));
rc = zeros( diff(v(1:2)), 1);

for k = 1:len_v
    rc = rc + cceps( ydrr(vidx) );    
end
rc = rc/len_v;

II = vidx;
plot(zca([ydry(II), ydrr(II) 0*icceps(cceps(ydrr(II))-rc)]))


