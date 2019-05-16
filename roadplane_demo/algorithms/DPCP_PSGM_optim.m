function [t, B, distances, i] = DPCP_PSGM_optim(Xtilde, c, mu_min, beta, maxiter)
t_start = tic;

if nargin < 3
    mu_min = 1e-15;
    maxiter = 200;
end

if nargin < 4
    maxiter = 200;
end

mu_0 = 1e-2;
alpha = 1e-3;

[D, N] = size(Xtilde);
obj = @(b)sum(abs(Xtilde'*b));

% initialization
% [B_0,~] = eigs(Xtilde*Xtilde',c,'SM');
[B_0, diag_0] = eig(Xtilde*Xtilde');
[~, ind] = sort(diag(diag_0));
B_0 = B_0(:, ind(1:c));
%bo = normc(randn(D,1));

for j = 1:c
    i = 1;
    b = B_0(:, j);
    mu = mu_0;
    if j == 1
        obj_old = obj(b);
        while mu > mu_min && i <= maxiter
            i = i + 1;
            %             grad = sum( sign(b'*Xtilde).*Xtilde, 2);
            grad = Xtilde * sign(Xtilde'*b);
            grad_norm = norm(grad)^2;
            %%% line search
            bk = b - mu * grad;
            while (obj(bk./norm(bk)) > obj_old - alpha * mu * grad_norm) && mu > mu_min;
                mu = mu * beta;
                bk = b - mu * grad;
            end
            b = bk ./ norm(bk);
            obj_old = obj(b);
        end
    else
        b = normc(b-B*(B' * b));
        obj_old = obj(b);
        while mu > mu_min && i <= maxiter
            i = i + 1;
            
            %%% line search
            %             grad = sum( sign(b'*Xtilde).*Xtilde, 2);
            grad = Xtilde * sign(Xtilde'*b);
            grad = grad - B * (B' * b);
            grad_norm = norm(grad)^2;
            bk = b - mu * grad;
            while (obj(bk./norm(bk)) > obj_old - alpha * mu * grad_norm) && mu > mu_min;
                mu = mu * beta;
                bk = b - mu * grad;
            end
            b = bk ./ norm(bk);
            obj_old = obj(b);
        end
    end
    B(:, j) = b;
    
end
t = toc(t_start);

distances = zeros(1, N);
for j = 1:N
    distances(j) = norm(B'*Xtilde(:, j));
end

end