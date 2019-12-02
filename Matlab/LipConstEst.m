function [L_est] = LipConstEst(x,A,maxit)
% Estimate Lipschitz constant

AT = @(x) A'*x;
A = @(x) A*x;
for j = 1:maxit
    z = A(x)/norm(x);
    x = AT(z);
end
L_est = norm(x);

end
