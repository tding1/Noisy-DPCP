function [inliers, inner, B] = distfn(B, X, t)

[D, N] = size(X);
inner = zeros(1,N);
inner = vecnorm(B' * X, 2, 1);
inliers = find(inner <= t);

end

