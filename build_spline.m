function [M, pp_spline] = build_spline(x, C_fd, h)
% Natural cubic spline through (x, C_fd)

    x   = x(:).';      % row
    Cfd = C_fd(:);     % column
    N   = length(x) - 1;

    M = zeros(N+1, 1);

    A_M = zeros(N-1, N-1);
    b_M = zeros(N-1, 1);

    for i = 1:N-1
        if i > 1
            A_M(i, i-1) = h;
        end
        A_M(i, i) = 4*h;
        if i < N-1
            A_M(i, i+1) = h;
        end
        b_M(i) = (6/h) * (Cfd(i+2) - 2*Cfd(i+1) + Cfd(i));
    end

    M(2:N) = gauss_naive(A_M, b_M);

    nIntervals = N;
    coefs = zeros(nIntervals,4);

    for i = 1:nIntervals
        a_i = (M(i+1)-M(i))/(6*h);
        b_i = M(i)/2;
        c_i = (Cfd(i+1) - Cfd(i))/h - (2*h*M(i) + h*M(i+1))/6;
        d_i = Cfd(i);
        coefs(i,:) = [a_i b_i c_i d_i];
    end

    pp_spline = mkpp(x, coefs);
end
