%=============================================
% For given a and B, solve
%     min_x  x'Bx - 2a'x  subject to ||x||_2=1, B>0

% =============================================

function [x_opt]=Gen_Quadratic_Constr_LS(a,B)
%normB = norm(a,'fro'); B = B/normB; a = a/normB;
%  Inputs --      a - N x 1        B - N x N
[N,~] = size(B);
[U,Sigma]=eig(B);
if max(abs(a)) < 1e-4
    [value,index] = min(diag(Sigma));
    x_opt = U(:,index);
else
    index = find(diag(Sigma)>1e-4); % ignore too small singular values and 0
    epsilon = (abs(U(:,index)'*a)).^2; % this is the epsilon^2 of the paper
    
    sigma = (diag(Sigma(index,index))); % this is the sigma^2 of the paper
    
    index = find(epsilon>1e-4);
    if size(index,1)==0
        [value,index] = min(diag(Sigma));
        x_opt = U(:,index);
        return
    end
    sigma = sigma(index); epsilon = epsilon(index);
    clear U Sigma index
    
    myfun = @(lambda) sum(epsilon./(lambda + sigma).^2)-1;
    lambda_up = sqrt(sum(epsilon)); % the upper bound
    lambda_lower = max(-sigma + sqrt(epsilon)); % the lower bound
    lambda_opt = fzero(myfun,[lambda_lower lambda_up]); % finding the largest root
    
    x_opt = pinv(lambda_opt*eye(N) + B)*a;
end


