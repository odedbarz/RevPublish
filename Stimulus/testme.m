function err = testme(x)

xa = hilbert(x);

xbar = real(abs(xa) .* exp(1j*angle(xa)));

err = norm(x - xbar);

plot([x, xbar, abs(xa)])
