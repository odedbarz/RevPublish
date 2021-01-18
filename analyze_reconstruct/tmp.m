
clc

n_bands = size(X_dry, 1);
n_lags = 10;

Adry = zca(im2col(X_dry, [n_bands, n_lags], 'sliding'));
Atst = zca(im2col(X_drr, [n_bands, n_lags], 'sliding'));
Aest = zca(im2col(X_est, [n_bands, n_lags], 'sliding'));



%%
% plot( diag(zca(Adry)'*zca(Aest)) )

a = Adry(:,nn);
b = Atst(:,nn);
c = Aest(:,nn);
cc_ac = corrcoef(a, c);     cc_ac = cc_ac(1,2);
cc_bc = corrcoef(b, c);     cc_bc = cc_bc(1,2);
cc_ab = corrcoef(b, a);     cc_ab = cc_ab(1,2);
cc_all_dry_est = corrcoef(X_dry(:), X_est(:));  cc_all_dry_est = cc_all_dry_est(1,2);
cc_all_drr_est= corrcoef(X_drr(:), X_est(:));   cc_all_drr_est = cc_all_drr_est(1,2);

fprintf('--> CC(dry-est): %g\n', cc_ac);
fprintf('--> CC(tst-est): %g\n', cc_bc);
fprintf('--> CC(dry-tst): %g\n', cc_ab);
fprintf('--> CC_ALL(dry-est): %g\n', cc_all_dry_est);
fprintf('--> CC_ALL(tst-est): %g\n', cc_all_drr_est);


fprintf('\n =========================\n');
nn = 73                        % worse correlation
nn = 253 %20 %90 %150 %420      % better correlation
Rx = @(x) reshape(x, [n_bands, n_lags]);
imagesc([Rx(Adry(:,nn)), Rx(Aest(:,nn)), Rx(Atst(:,nn))]);
title(  sprintf('$cc_(dry/estc)$: %.2f $----$ $cc_(drr/est)$: %.2f', cc_ac, cc_bc) );








