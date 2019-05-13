function [b, f] = DPCP_PSGM(X,O, noise)

% solves min_b ||Xtilde^T b||_1 s.t. b^T b=1
% INPUT
% Xtilde  : DxN data matrix of N data points of dimension D
% c        : Dimension of the orthogonal complement of the subspace.
%            To fit a hyperplane use c=1.
% mu_min  : minimum step size. Typicall is set to 10^(-15)
% maxiter : Maximal number of iterations. Typically set to 200.

% Parameters:
%alpha and beta: parameters for linear search. Typically is set to
%                  alpha = 0.001 and beta = 1/2.
% mu_0    : Initialization. Typically is set to 10^(-2).
%%%%
% OUTPUT
% distances: Distance of each point to the estimated subspace.
% b        : Dx1 vector in  the orthogonal complement of the subspace.

if nargin < 3
   noise = 0; 
end

mu_min = 1e-15; maxiter = 400;

Xtilde = [X + noise, O];
% Xtilde = normc(Xtilde);

mu_0 = 1e-2; alpha = 1e-3; beta = 1/2;

obj = @(b)norm(Xtilde'*b,1);

% initialization
[b,~] = eigs(Xtilde*Xtilde',1,'SM');

%bo = normc(randn(D,1));
mu = mu_0;
obj_old = obj(b);
i = 1;
while mu>mu_min && i<= maxiter
    i = i+1;
    
%     grad = sum( repmat(sign(b'*Xtilde),D,1).*Xtilde, 2);
    grad = Xtilde*sign(Xtilde'*b);
    grad_norm = norm(grad)^2;
    % line search
    tmp = b - mu*grad;
    while (obj( tmp / norm(tmp) )> obj_old - alpha*mu*grad_norm) && mu>mu_min
        mu = mu*beta;
        tmp = b - mu*grad;
    end
    b = tmp / norm(tmp);
    obj_old = obj(b);
end

f = obj(b);

end















