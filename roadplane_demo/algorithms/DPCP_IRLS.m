function [B, t_elapsed, k] = DPCP_IRLS(X_tilde, c, delta, T, epsilon_J, budget)
t_start = tic;
[D, N] = size(X_tilde);
Delta_J = Inf;
k = 0;
w = ones(1, N);
J_old = inf(1, N);
J_new = zeros(1, N);

while (Delta_J > epsilon_J) && (k < T) && toc(t_start) <= budget
    dw = sparse(1:N, 1:N, w);
    R_X = X_tilde * dw * X_tilde';
    %     [U, ~, ~] = svd(R_X);
    %     B = U(:,D-c+1:D);
    [U, L] = eig(R_X);
    lambda = diag(L);
    [lambda, perm] = sort(lambda, 'descend');
    B = U(:, perm(end));
    J_new = vecnorm(B'*X_tilde, 2, 1);
    w = 1 ./ max(J_new, delta);
    k = k + 1;
    Delta_J = 1 - sum(J_new) / (sum(J_old) + 10^(-9));
    J_old = J_new;
end

t_elapsed = toc(t_start);

% disp(['DPCP-IRLS finished with delta_J: ', num2str(Delta_J)]);
% disp(['DPCP-IRLS iter: ', num2str(k), ' time: ', num2str(t_elapsed), ' t/k: ', num2str(t_elapsed/k)]);


end
