function TV2_rec = Solve_TV2(sino, A, N, alpha)


f = zeros(N*N, 1);
lambda = 10^(-5);
beta = 10^-5;
maxit = 200;
sino = sino(:);

grad_diff = 2*A.'*(A*f(:) - sino);
gradasd = gradjubs(f, beta, N);
grad_funct_now = grad_diff + alpha*gradasd;
iter = 1;
while iter <= maxit
    fnext = f - lambda * grad_funct_now;
    fnext(fnext<0) = 0;
    grad_diff = 2*A.'*(A*fnext(:) - sino);
    gradasd = gradjubs(fnext, beta, N);
    grad_funct_next = grad_diff + alpha*gradasd;
    
    if iter < maxit
        lambda = ((fnext - f).'*(fnext - f)) / ((fnext - f).'*(grad_funct_next - grad_funct_now));
    end
    
    f = fnext;
    grad_funct_now = grad_funct_next;
    iter = iter + 1;
    disp(iter);
end

TV2_rec = reshape(f, [N,N]);
end