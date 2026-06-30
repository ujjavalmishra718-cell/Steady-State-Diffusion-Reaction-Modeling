function plot_results(x, C_fd, C_exact, pp_spline, L)

    x   = x(:).';
    Cfd = C_fd(:);
    Cex = C_exact(:);

    xf = linspace(0, L, 400);
   
    C_spline_fine = ppval(pp_spline, xf);

    figure;
    hold on; grid on;

    plot(x, Cfd, 'r', 'MarkerSize', 4, 'DisplayName', 'Numerical (FD)');
    plot(x, Cex, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
   
    plot(xf, C_spline_fine, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Spline');

    xlabel('x (cm)');
    ylabel('C(x)');
    title('Comparison: FD, Spline, and Analytical Solution');
    legend('Location', 'best');

    error_fd = Cfd - Cex;
    fprintf('   i      x_i      C_fd        C_exact      error (FD - exact)\n');
    N = length(x) - 1;
    for i = 1:N+1
        fprintf('%4d   %6.3f   %9.5f   %9.5f   %+10.5e\n', ...
                i-1, x(i), Cfd(i), Cex(i), error_fd(i));
    end
end
