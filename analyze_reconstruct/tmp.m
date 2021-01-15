
clc

n_bands = size(X_test0, 1);
n_lags = 10;

X_est = obj.X_est;

Adry = zca(im2col(X_test0, [n_bands, n_lags], 'sliding'));
Atst = zca(im2col(X_test_kn, [n_bands, n_lags], 'sliding'));
Aest = zca(im2col(X_est, [n_bands, n_lags], 'sliding'));



%%
plot( diag(zca(Adry)'*zca(Aest)) )

fprintf('\n =========================\n');
nn = 310                        % worse correlation
nn = 310 %20 %90 %150 %420      % better correlation
Rx = @(x) reshape(x, [n_bands, n_lags]);

imagesc([Rx(Adry(:,nn)), Rx(Aest(:,nn)), Rx(Atst(:,nn))]);

a = Adry(:,nn);
b = Atst(:,nn);
c = Aest(:,nn);
cc_ac = corrcoef(a, c);
cc_bc = corrcoef(b, c);
cc_ab = corrcoef(b, a);
fprintf('--> CC(dry-est): %g\n', cc_ac(1,2));
fprintf('--> CC(tst-est): %g\n', cc_bc(1,2));
fprintf('--> CC(dry-tst): %g\n', cc_ab(1,2));


corrcoef(X_test0(:), X_est(:))
corrcoef(X_test_kn(:), X_est(:))






