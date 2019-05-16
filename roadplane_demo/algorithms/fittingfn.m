function B = fittingfn(X,dim)

[U, S, V] = svd(normc(X));
B = U(:,dim+1:end);

end

