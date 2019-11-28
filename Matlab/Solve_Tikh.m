function [Tik_rec] = Solve_Tikh(sino, A)
    maxit = 200;
    f0 = zeros(512);
    lambda_min = 1e-10;
    lambda_max = 1e10;
    gamma = 1.e-4;
    beta = 0.4;
    lambda = 1.3;
    alpha = 1;
    stopcrit = 1;
    tolstop = 1e-4;

    f_size = size(f0);
    f = f0(:);
    
    sino = sino(:);
    
    funval = 0.5 * ((norm(A*f-sino))^2 + alpha * (norm(f))^2);
    g = A' * (A*f-sino) + alpha * f;
    funct = zeros(maxit + 1, 1);
    funct(1) = funval;
        
    iteration = 1;
    notDone = true;
    
    while notDone
        fprintf('Iteration %d\n', iteration);

        z = f - lambda * g;
        z(z < 0) = 0;

        d = z - f;
        gd = dot(d, g);

        mu = 1;

        fcontinue = 1;

        fr = funval;

        while fcontinue
            fplus = f + mu*d;

            funval = 0.5 * ( (norm(A*fplus-sino))^2 + alpha * (norm(fplus))^2 );

            % Sufficient decrease condition
            if ( funval <= fr + gamma * mu * gd || mu<1e-12)
                f = fplus; clear fplus;
                sk = mu*d;
                gtemp = A' * (A*f-sino) + alpha * f;
                zetak = gtemp - g;
                g = gtemp; clear gtemp; 
                fcontinue = 0;
            else
                mu = mu * beta;
            end
        end

        if (funval >= fr)
            disp('Warning: function value >= fr');
        end

        bk = dot(sk(:),zetak(:)); 

        % BB1 updating rule 
        lambda = dot(sk(:),sk(:))/bk;

        % BB2 updating rule
        % lambda = bk/ dot(zetak(:),zetak(:));

        lambda = min( max(lambda, lambda_min), lambda_max );
        % lambda = double(single(lambda));

        funct(iteration+1) = funval;

        % Stop criteria
        switch stopcrit
            case 1
                % Maximum number of iterations
            case 2
                % || x_k - x_k-1 ||_2 / || x_k ||_2
                crit_value = norm(sk)/norm(f);
                notDone = crit_value > tolstop;
            case 3
                % | f_k - f_k-1 | / |f_k|
                crit_value = abs(funval - funct(iteration))/abs(funval); 
                notDone = crit_value > tolstop;
         end

        if iteration >= maxit
                notDone = false;
        end       
       
    iteration = iteration + 1;
    end
    
    Tik_rec = reshape(f, f_size);
    iteration = iteration - 1;
    functs = funct(1:iteration);
    figure;
    imshow(Tik_rec, [])
end