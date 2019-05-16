function [result] = metric(B_hat, X_tilde, L, B_gt, th, points3D, vis)
P = abs(B_hat'*X_tilde);

if vis == 1
    [AUCPR, AUCROC] = prMetricArea(P, L, 1);
    n_points = length(X_tilde);
    inliers_hat = find(P < th);
    outliers_hat = setdiff(1:n_points, inliers_hat);
    figure;
    title('Estimated, using ground-truth threshold');
    pcshow(points3D(:, inliers_hat)', [1, 0, 0]);
    hold on;
    pcshow(points3D(:, outliers_hat)', [0, 0, 1]);
    hold off;
else
    [AUCPR, AUCROC] = prMetricArea(P, L, 0);
end

B_gt_denormalized = B_gt ./ norm(B_gt(1:3));
B3_gt = B_gt_denormalized(1:3);
gamma_gt = B_gt_denormalized(4);

B_hat_denormalized = B_hat ./ norm(B_hat(1:3));
B3_hat = B_hat_denormalized(1:3);
gamma_hat = B_hat_denormalized(4);


theta4 = acos(abs(B_hat'*B_gt))/pi*180;
theta3 = acos(abs(B3_hat'*B3_gt))/pi*180;
error_t = abs((-gamma_gt+B3_gt'*B3_hat*gamma_hat)/gamma_gt);

result.AUCPR = AUCPR;
result.AUCROC = AUCROC;
result.theta4 = theta4;
result.theta3 = theta3;
result.B = B_hat;
result.error_t = error_t;
end