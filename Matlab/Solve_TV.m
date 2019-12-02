function [TV_rec] = Solve_TV(sino, A, N, alpha)

% Input: A = system matrix
%        sino = sinogram
%        W = matrix of weights
%        x0 = initial guess
%        obj = true solution
%
% Output: x = solution
%         iter = number of iteration to convergence
%         err = vector containing the error at each iteration
%         funct = vector containing the objective function value at each iteration 


maxit = 100;
maxpowit = 10;
x = zeros(N);
x = x(:);
gamma = 1.e-4;
stopcrit = 3;
tolstop = 1e-2;

% Allocation
iter = 1;
loop = true;
err = [];
funct = [];

% Vectorization
sino = sino(:);

D = spdiags([-ones(N-1,1), ones(N-1,1); 0, 1], 0:1, N, N);
D1 = kron(speye(N),D);
D2 = kron(D, speye(N));
Dhuge = [D1; D2];
D1 = @(x) D1*x;
D2 = @(x) D2*x;
regu = norm(sqrt(D1(x).^2 + D2(x).^2),1);

matrix = [gamma*Dhuge; A];
ell0 = randn(size(matrix,2),1);
[L_est] = LipConstEst(ell0,matrix,maxpowit);
L = normest(matrix)^2;
fprintf('\n L = %f \t L_est=%f ', L, L_est);

AT = @(x) A'*x; 
A = @(x) A*x;

tau = alpha*0.99/(sqrt(L_est));
sigma = 0.99/(alpha*sqrt(L_est));

y = zeros(2*N^2,1); % initial y
z = zeros(size(sino)); % initial z

%%%%%%%%%%%%%%%%%%

LS = A(x) - sino;
fv = 0.5*(norm(LS)^2 + gamma*regu);
funct(1) = fv;

while loop
    % step x
    tmp = gamma*Dhuge'*y + AT(z);
    x_hat = x - tau*tmp;
    x_hat(x_hat < 0) = 0;
    x_new = x_hat;
    
    % step y
    proxgs = @(y) y./repmat(max(1,sqrt(sum(reshape(y,N^2,2).^2,2))),2,1);
    y_new = proxgs(y + sigma*gamma*Dhuge*(2*x_new-x));
    
    % step z
    z_hat = z + sigma*A(2*x_new-x);
    z_new = (z_hat - sigma*sino)/(1+sigma);
    
    
    % stopping criterion    
    LStry = A(x_new)-sino;
    ftry = 0.5*(norm(sqrt(LStry))^2 + gamma*(norm(sqrt(D1(x_new).^2 +D2(x_new).^2),1)));
    
    if iter >= 2
        switch stopcrit
            case 1
                normstep_x = norm(x_new-x)/norm(x);
                normstep_y = norm(y_new-y)/norm(y);
                normstep_z = norm(z_new-z)/norm(z);
                normstep = normstep_x + normstep_y + normstep_z;
                fprintf('\n \t || x_k - x_(k-1) ||^2 / || x_k ||^2 %e tol %e\n',normstep,tolstop);
                loop = ( normstep > tolstop );
            case 2
                prim = - gamma*Dhuge'*y_new - AT(z_new);
                prim( x==0 ) = max(0,prim(x==0)); 
                
                Y = reshape(y_new,[],2);
                DX = reshape(Dhuge*x_new,[],2);
                normy = sqrt(sum(Y.^2,2));
                normDX = sqrt(sum(DX.^2,2));
                tmp = (abs(1-normy) <= 1e-6).*normDX;
                dual = Y.*repmat(tmp,1,2)-DX;
                dual = sqrt(sum(dual.^2,2)); % vector of norms
                
                pd = A(x_new)-sino-z_new;
                
                normstep = max(norm(prim),max(norm(dual),norm(pd)));
                fprintf('\n \t primal-dual %e tol %e\n',normstep,tolstop);
                loop = ( normstep > tolstop );
        end
    end
    iter = iter + 1;
    x = x_new;
    y = y_new;
    z = z_new;
    
    funct(iter) = ftry;
    
    if iter >= maxit
        loop = false;
    end 
end
iter = iter - 1;
funct = funct(1:iter);

TV_rec = x;

end