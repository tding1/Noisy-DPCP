function [resultRANSAC] = runRansac(X_tilde, th_gt, B_gt, L, budget, nExp, points3D, vis)
thArray = 10.^(-3:0.1:0);
AUCArray = [];
kArray = [];
tArray = [];
theta4Array = [];
theta3Array = [];
error_tArray = [];
f1sArray = [];
for iExp = 1:nExp
    tprArray = [];
    fprArray = [];
    for iTh = 1:length(thArray)
        th = thArray(iTh);
        [B_RANSAC, ~, best_inliers, kArray(end+1), tArray(end+1)] = RANSAC_PCA(X_tilde, 3, th, budget);
        [tprArray(end+1), fprArray(end+1)] = tprFpr(L, best_inliers);
    end
    [fprArray, indices] = sort(fprArray);
    tprArray = tprArray(indices);
    AUCArray(end+1) = trapz([0, fprArray, 1], [0, tprArray, 1]);
    [B_RANSAC, ~, ~, ~, ~] = RANSAC_PCA(X_tilde, 3, th_gt, budget);
    [resulttemp] = metric(B_RANSAC, X_tilde, L, B_gt, th_gt, points3D, vis);
    theta4Array(end+1) = resulttemp.theta4;
    theta3Array(end+1) = resulttemp.theta3;
    error_tArray(end+1) = resulttemp.error_t;
end
resultRANSAC = struct();
resultRANSAC.AUCPR = nan;
resultRANSAC.AUCROC = mean(AUCArray);
resultRANSAC.theta4 = mean(theta4Array);
resultRANSAC.theta3 = mean(theta3Array);
resultRANSAC.error_t = mean(error_tArray, 2);
resultRANSAC.B = B_RANSAC;
resultRANSAC.k = mean(kArray);
resultRANSAC.t = mean(tArray);
end
