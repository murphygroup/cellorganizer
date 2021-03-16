function [X, f_mat] = box_constrained_lsq_admm(A, b, ub, lb, nIter, rho)

if nargin < 5
    nIter = 15;
end
if nargin < 6
    rho = 1e-6;
end

[n, p] = size(A);

% function handles for the objective function and prox operator
f_func = @(x) norm(A * x - b);
prox_func = @(x, t) max(0, x) + t; 

X = rand(p, 1);
U_1 = zeros(p, 1);
U_2 = zeros(p, 1);
Z_1 = X;
Z_2 = -X;

% gradient function handle. 

f_mat = zeros(nIter + 1, 1);
f_mat(1) = f_func(X);
% methods. 
ATA_reg = (A' * A + 2 * rho * speye(p));
ATb = A' * b;
for k = 1 : nIter

    X = ATA_reg \ (ATb + rho * ( Z_1 + U_1 - Z_2 - U_2));
    Z_1 = prox_func(X - U_1 - lb, lb);
    Z_2 = prox_func(-X - U_1 + ub, -ub);

    U_1 = U_1 + Z_1 - X;
    U_2 = U_2 + Z_2 + X;

    f_k = f_func(X);
    f_mat(k + 1) = f_k;
end


end


