function [x, w] = gauss_legendre_1d(n)
%GAUSS_LEGENDRE_1D Gauss-Legendre points and weights on [-1, 1].

    if n == 1
        x = 0;
        w = 2;
        return;
    end

    k = (1:n-1).';

    beta = k ./ sqrt(4*k.^2 - 1);

    J = diag(beta, 1) + diag(beta, -1);

    [V, D] = eig(J);

    x = diag(D);

    [x, order] = sort(x);

    w = 2 * (V(1, order).^2);

    x = x(:);
    w = w(:);
end