function gradasd = gradjubs(f, beta, n)
N = n*n;
gradasd = zeros(N, 1);
f = f(:);
for iii = 1:N
    if iii + n <= N
        fplusn = f(iii + n);
    else
        fplusn = 0;
    end
    if iii - n >= 1
        fminusn = f(iii - n);
    else
        fminusn = 0;
    end
    if iii + n - 1 <= N
        fplusn1 = f(iii + n - 1);
    else
        fplusn1 = 0;
    end
    if iii - n + 1 >= 1
        fminusn1 = f(iii - n + 1);
    else
        fminusn1 = 0;
    end
    if iii + 1 <= N
        fplus1 = f(iii + 1);
    else
        fplus1 = 0;
    end
    if iii - 1 >= 1
        fminus1 = f(iii - 1);
    else
        fminus1 = 0;
    end
    
    gradasd(iii) = ((2*f(iii) - fplusn - fplus1) / sqrt((fplusn - f(iii))^2 + (fplus1 - f(iii))^2 + beta)) + ((f(iii) - fminusn) / sqrt((f(iii) - fminusn)^2 + (fminusn1 - fminusn)^2 + beta)) + ((f(iii) - fminus1) / sqrt((f(iii) - fminus1)^2 + (fplusn1 - fminus1)^2 + beta));
end
end