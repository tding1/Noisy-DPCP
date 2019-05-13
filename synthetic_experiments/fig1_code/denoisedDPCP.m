function [b, y, f] = denoisedDPCP(X,O, noise, tau)

maxiter = 400;
Xtilde = [X + noise, O];
% Xtilde = normc(Xtilde);

% [U,D,~] = svd(Xtilde * Xtilde');

% sz = size(X, 1);
% inv_xxt = Xtilde * Xtilde' \ eye(sz);
[b,~] = eigs(Xtilde*Xtilde',1,'SM');
tmp = Xtilde' * b;
y = max(abs(tmp)-tau, 0) .* sign(tmp);

obj = @(b, y) tau * norm(y,1) + 0.5 * norm(Xtilde'*b - y)^2;
old_obj = obj(b, y);
fval = zeros(maxiter,1);
fval(1) = old_obj;

delta_obj = inf;
epi = 1e-6;

i = 1;
while i <= maxiter && delta_obj > epi * old_obj
    tmp = Xtilde' * b;
    y = max(abs(tmp)-tau, 0) .* sign(tmp);
    b = Gen_Quadratic_Constr_LS(Xtilde*y, Xtilde*Xtilde');

    i = i + 1;
    new_obj = obj(b,y);
    delta_obj = old_obj - new_obj;
    old_obj = new_obj;
    fval(i) = new_obj;
    
%     fbb = obj(bb, y);
end

% plot(fval)

f = obj(b,y);


end

