function [tL21, e, normal, iter] = RPCA_L21_OutlierDetection(Y, tau, lambda, budget)

tstart = tic;
% L21 robust PCA
[L, E, iter] = rpca(Y, 'L21', tau, lambda, budget);
% disp(['L21 runs ', num2str(iter), ' iterations']);

[U, ~, ~] = svd(L*L');
normal = U(:, end);

e = Y' * normal;

tL21 = toc(tstart);

end
