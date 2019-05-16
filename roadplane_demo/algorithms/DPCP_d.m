function [y, b, counter] = DPCP_d(X_tilde,n_0,T_ADM,epsilon_J,tau,local_budget)

t_start = tic;
% initializations
[D, N] = size(X_tilde);
DJ = inf;
b = n_0;
y = X_tilde'*b;
J_old = tau*norm(y,1)+ 0.5*norm(y-X_tilde'*b,2)^2;
counter = 0;
RX_inv_X = (X_tilde*X_tilde'+10^(-6)*eye(D)) \ X_tilde;
while (counter < T_ADM) && (DJ > epsilon_J) && (toc(t_start)<local_budget)
    counter = counter + 1;
    
    % update y            
    y = wthresh(X_tilde'*b,'s',tau);
    
    % update b
    if norm(y) == 0
        break;
    else
        b = RX_inv_X *y;
        b = b / (norm(b)+10^(-6));
    end

    
    % updates
    J_new = tau*norm(y,1) + 0.5*norm(y-X_tilde'*b,2)^2;
    DJ = (J_old - J_new) / (J_old+10^(-6));         
    J_old = J_new;                                 
end
timeAlgo = toc(t_start);
% disp(['DPCP-d iter: ', num2str(counter), ' time: ', num2str(timeAlgo), ' t/k: ', num2str(timeAlgo/counter)]);

end

