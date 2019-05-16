function P = REAPER_IRLS_step(X,w,c)

[D,N] = size(X);
nu = zeros(D,1);
d = D - c;

dw = sparse(1:N, 1:N, w);
R_X = X * dw * X';
[U, L] = eig(R_X);
lambda = diag(L);
[lambda perm] = sort(lambda,'descend');
U = U(:,perm);

if lambda(d+1) == 0
    nu(1:d) = ones(d,1);
else    
    for i = d+1 : D
        theta = (i-d) / sum(ones(i,1) ./ lambda(1:i));
        if (i<D) && (lambda(i)>theta) && (theta>=lambda(i+1))
            break;
        end
    end    
    for i = 1 : D    
        if lambda(i)>theta
            nu(i) = 1 - (theta/lambda(i));
        else 
            nu(i) = 0;
        end
    end
end

P = U * diag(nu) * U';

end

