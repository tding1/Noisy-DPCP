function [t, distances, B, k] = REAPER_IRLS_optim(X,c,delta,T_max,epsilon_J,budget)

t_start = tic;
[D, N] = size(X);
Delta_J = Inf;
k = 0;
w = ones(1, N);
J_old = zeros(1, N);
J_new = zeros(1, N);
while (k<T_max) && (Delta_J>epsilon_J) && (toc(t_start) < budget)
    P = REAPER_IRLS_step_optim(X,w,c);
    J_new = sqrt(sum((X-P*X).^2, 1));
    w = 1./max(J_new, delta);
    k = k + 1;
    Delta_J = abs(sum(J_old)-sum(J_new))/(sum(J_old)+10^(-9));
    J_old = J_new;
end

[U S V] = svd(P);
B = U(:,D-c+1:D);
t = toc(t_start);
% disp(['REAPER finished in ', num2str(k), ' iterations'])

if c == 1
    distances = abs(B'*X);
else
    for j = 1 : N
        distances(j) = norm(B'*X(:,j));
    end
end

end

