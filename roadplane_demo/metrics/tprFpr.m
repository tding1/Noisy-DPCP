function [tpr, fpr] = tprFpr(L, identified)
P = find(L == 1);
TP = intersect(P, identified);
tpr = length(TP) / length(P);
N = find(L == 0);
FP = intersect(N, identified);
fpr = length(FP) / length(N);

end

