clear; clc; close all;

%% PARAMETERS
L  = 5.0;        % cm
h  = 0.05;       % grid spacing, cm
k  = 0.1;        % day^-1
D0 = 0.84;       % cm^2/day at x=0
DL = 1.09;       % cm^2/day at x=L

m  = (DL - D0)/L;        % slope of D(x)  => here 0.05
mh_over2 = m*h/2;        % m*h/2 (used in FD coefficients)
hk       = h*k;

% grid
x = 0:h:L;               % x_0 ... x_N
N = length(x) - 1;       % N = 100, last index corresponds to x = L

%% ANALYTICAL SOLUTION
F       = log(DL/D0);
A_const = -100 * k / (1 + (k/m)*F);
C_exact = 100 + (A_const/m) * log( (D0 + m*x) / D0 );

%% FINITE DIFFERENCE SYSTEM

C0 = 100;
Afd = zeros(N+1, N+1);
b   = zeros(N+1, 1);

% i = 0: C_0 = 100
Afd(1,1) = 1;
b(1)     = C0;

% interior nodes: i = 1..N-1
for i = 1:N-1
    xi = x(i+1);
    Di = D0 + m*xi;

    left  = Di - mh_over2;
    cent  = -2*Di;
    right = Di + mh_over2;

    row = i+1;
    Afd(row, i  ) = left;
    Afd(row, i+1) = cent;
    Afd(row, i+2) = right;
end

% last node: i = N
xN = x(end);
DN = D0 + m*xN;
coef_CNm1 = DN;
coef_CN   = -DN - hk*(1 + (m*h)/(2*DN));

rowN = N+1;
Afd(rowN, N  ) = coef_CNm1;
Afd(rowN, N+1) = coef_CN;

% solve with naive Gauss
C_fd = gauss_naive(Afd, b);    % [C_0; ...; C_N]

%% NATURAL CUBIC SPLINE
[M, pp_spline] = build_spline(x, C_fd, h);

%% PLOT COMPARISONS + ERROR TABLE
plot_results(x, C_fd, C_exact, pp_spline, L);
