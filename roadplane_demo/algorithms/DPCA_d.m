function [B, timeAlgo, k_all] = DPCA_d(X_tilde,T_ADM,epsilon_J,c,tau,NORMALIZATION,global_budget)

t_start = tic;
% dimensions/initializations
[D, N] = size(X_tilde);
B = zeros(D,c);
Pi = eye(D);
Y = zeros(N,c);
k_all = 0;
% find c orthogonal dual principal components
for i = 1 : c
    Xp = Pi'*X_tilde;
    if NORMALIZATION == 1
        Xp = normc(Xp);
    end
    % find a dpc
    [V, d] = eig(Xp*Xp');
    [d,ind] = sort(diag(d), 'descend');
    V = V(:, ind);
    
    temp = (global_budget-toc(t_start))/(c-i+1);
    if temp > 0
        local_budget = temp;
    else
        local_budget = 0.001;
    end
    
    [y, b, k] = DPCP_d(Xp,V(:,end),T_ADM,epsilon_J,tau,local_budget);
    k_all = k_all + k;
    
    Y(:,i) = abs(Xp'*b);
    B(:,i) = Pi * b;
    % update the projection Pi
    [Q R] = qr([B(:,1:i) eye(D)]);
    Pi = Q(:,i+1:end);
end

timeAlgo = toc(t_start);
% disp(['DPCP-d finished in ', num2str(timeAlgo), 's']);
end

