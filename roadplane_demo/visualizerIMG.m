function [] = visualizerIMG(date, seq, frame, dataPath, resultsAll)
W = 1243;
H = 376;
n_row = 5;
n_col = 2;
fontColor = [1, 1, 1];
fontsize = 20;
left_b = 10;
top_b = 330;
counter = 1;

imgPath = [dataPath, sprintf('09_%02d_%03d/%02d_%03d_%04d.png', date, seq, date, seq, frame)];
img = imread(imgPath);

load([dataPath, sprintf('09_%02d_%03d/%02d_%03d_%04d_data.mat', date, seq, date, seq, frame)]);


[U, ~, ~] = svd(normc(inliers), 'econ');
B_GT = U(:, end);

gt_threshold = max(abs(B_GT'*normc(inliers)));

B_D = resultsAll.D.B;
B_IRLS = resultsAll.IRLS.B;
B_L21 = resultsAll.L21.B;
B_PSGM = resultsAll.PSGM.B;
B_REAPER = resultsAll.REAPER.B;
B_SVD = resultsAll.SVD.B;
B_RANSAC = resultsAll.RANSAC.B;
B_RANSAC_10 = resultsAll.RANSAC_10.B;
B_RANSAC_100 = resultsAll.RANSAC_100.B;

n_inliers = resultsAll.n_inliers;
n_outliers = resultsAll.n_outliers;

X = [inliers, outliers];
X_tilde = normc([inliers, outliers]);
target_gt = [ones(1, n_inliers), zeros(1, n_outliers)];
idI_gt = find(target_gt > 0.5);
idO_gt = find(target_gt < 0.5);

%% calculate the label for each method
[idI_D, idO_D] = cut(X_tilde, B_D, gt_threshold);
[idI_IRLS, idO_IRLS] = cut(X_tilde, B_IRLS, gt_threshold);
[idI_L21, idO_L21] = cut(X_tilde, B_L21, gt_threshold);
[idI_PSGM, idO_PSGM] = cut(X_tilde, B_PSGM, gt_threshold);
[idI_REAPER, idO_REAPER] = cut(X_tilde, B_REAPER, gt_threshold);
[idI_SVD, idO_SVD] = cut(X_tilde, B_SVD, gt_threshold);
[idI_RANSAC, idO_RANSAC] = cut(X_tilde, B_RANSAC, gt_threshold);
[idI_RANSAC_10, idO_RANSAC_10] = cut(X_tilde, B_RANSAC_10, gt_threshold);
[idI_RANSAC_100, idO_RANSAC_100] = cut(X_tilde, B_RANSAC_100, gt_threshold);

[im_GT, normal_GT] = (Draw_image(img, inliers, outliers, date, B_GT));
[im_D, normal_D] = (Draw_image(img, X(:, idI_D), X(:, idO_D), date, B_D));
[im_IRLS, normal_IRLS] = (Draw_image(img, X(:, idI_IRLS), X(:, idO_IRLS), date, B_IRLS));
[im_L21, normal_L21] = (Draw_image(img, X(:, idI_L21), X(:, idO_L21), date, B_L21));
[im_PSGM, normal_PSGM] = (Draw_image(img, X(:, idI_PSGM), X(:, idO_PSGM), date, B_PSGM));
[im_REAPER, normal_REAPER] = (Draw_image(img, X(:, idI_REAPER), X(:, idO_REAPER), date, B_REAPER));
[im_SVD, normal_SVD] = (Draw_image(img, X(:, idI_SVD), X(:, idO_SVD), date, B_SVD));
[im_RANSAC, normal_RANSAC] = (Draw_image(img, X(:, idI_RANSAC), X(:, idO_RANSAC), date, B_RANSAC));
[im_RANSAC_10, normal_RANSAC10] = (Draw_image(img, X(:, idI_RANSAC_10), X(:, idO_RANSAC_10), date, B_RANSAC_10));
[im_RANSAC_100, normal_RANSAC100] = (Draw_image(img, X(:, idI_RANSAC_100), X(:, idO_RANSAC_100), date, B_RANSAC_100));

im_GT = imresize(im_GT, [H, W]);
im_IRLS = imresize(im_IRLS, [H, W]);
im_D = imresize(im_D, [H, W]);
im_REAPER = imresize(im_REAPER, [H, W]);
im_L21 = imresize(im_L21, [H, W]);
im_SVD = imresize(im_SVD, [H, W]);
im_PSGM = imresize(im_PSGM, [H, W]);
im_RANSAC = imresize(im_RANSAC, [H, W]);
im_RANSAC_10 = imresize(im_RANSAC_10, [H, W]);
im_RANSAC_100 = imresize(im_RANSAC_100, [H, W]);

im_all = imresize(im_GT, [H * n_row, W * n_col]);
im_all(1:H, 1:W, :) = im_GT;
im_all(1:H, W+1:2*W, :) = im_PSGM;
im_all(H+1:2*H, 1:W, :) = im_IRLS;
im_all(H+1:2*H, W+1:2*W, :) = im_D;
im_all(2*H+1:3*H, 1:W, :) = im_SVD;
im_all(2*H+1:3*H, W+1:2*W, :) = im_REAPER;
im_all(3*H+1:4*H, 1:W, :) = im_L21;
im_all(3*H+1:4*H, W+1:2*W, :) = im_RANSAC;
im_all(4*H+1:5*H, 1:W, :) = im_RANSAC_10;
im_all(4*H+1:5*H, W+1:2*W, :) = im_RANSAC_100;
figure;
imshow(im_all);

text(left_b, top_b, 'GroundTruth', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b+W, top_b, 'DPCP-PSGM', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b, top_b+H, 'DPCP-IRLS', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b+W, top_b+H, 'DPCP-d', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b, top_b+H*2, 'SVD', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b+W, top_b+H*2, 'REAPER', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b, top_b+H*3, '$\ell_{2,1}$-RPCA', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b+W, top_b+H*3, 'RANSACx1', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b, top_b+H*4, 'RANSACx10', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
text(left_b+W, top_b+H*4, 'RANSACx100', 'FontSize', fontsize, 'Color', fontColor, 'Interpreter', 'latex');
counter = counter + 1;
end

function [Roadimage, normal_vec] = Draw_image(img, inliers, outliers, date, B)

%% Camera matrixs
center = [8.174; -0.11; -1.631; 1];
B = B ./ B(4);
end_normal = center + B;
end_normal = end_normal / end_normal(4);
if date == 29
    velo2rect = [0.0048523, -0.00302069, 0.99998367, 0.; ...
        -0.99992978, 0.01079808, 0.00488465, 0.; ...
        -0.01081266, -0.99993711, -0.00296808, 0.; ...
        -0.00711321, -0.06176636, -0.2673906, 1.];
    
    rect2disp = [718.3351, 0., 0., 0.; ...
        0., 718.3351, 0., 0.; ...
        600.3891, 181.5122, 0., 1.; ...
        0., 0., 385.8846, 0.];
else
    velo2rect = [2.34773698e-04, 1.04494074e-02, 9.99945389e-01, 0.00000000e+00; ...
        -9.99944155e-01, 1.05653536e-02, 1.24365378e-04, 0.00000000e+00; ...
        -1.05634778e-02, -9.99889574e-01, 1.04513030e-02, 0.00000000e+00; ...
        -2.79681694e-03, -7.51087914e-02, -2.72132796e-01, 1.00000000e+00];
    
    rect2disp = [721.5377, 0., 0., 0.; ...
        0., 721.5377, 0., 0.; ...
        609.5593, 172.854, 0., 1.; ...
        0., 0., 387.5744, 0.];
end

%% normals
center_img = center' * velo2rect * rect2disp;
center_img = center_img';
center_img = center_img / center_img(4);
end_normal_img = end_normal' * velo2rect * rect2disp;
end_normal_img = end_normal_img';
end_normal_img = end_normal_img / end_normal_img(4);

center2d_img = round(center_img(1:2)) + 1;
end_normal_img = round(end_normal_img(1:2)) + 1;

normal_vec = [center2d_img; end_normal_img];

%% inliers
inliers_img = inliers' * velo2rect * rect2disp;
inliers_img = inliers_img';
inliers_img = inliers_img ./ inliers_img(4, :);

inliers_visiable_idx = find(inliers_img(1, :) <= 1242-1 & ...
    inliers_img(1, :) >= 0 & ...
    inliers_img(2, :) <= 375-1 & ...
    inliers_img(2, :) >= 0 & ...
    inliers_img(3, :) <= 255 & ...
    inliers_img(3, :) >= 0);
inliers_img = inliers_img(:, inliers_visiable_idx);

inliers2d_img = round(inliers_img(1:2, :)) + 1;
for i = 1:size(inliers2d_img, 2)
    y = inliers2d_img(1, i);
    x = inliers2d_img(2, i);
    img(x, y, 1) = 0;
    img(x, y, 2) = 0;
    img(x, y, 3) = 255;
    
    img(x+1, y, 1) = 0;
    img(x+1, y, 2) = 0;
    img(x+1, y, 3) = 255;
    
    img(x, y+1, 1) = 0;
    img(x, y+1, 2) = 0;
    img(x, y+1, 3) = 255;
    
    img(x+1, y+1, 1) = 0;
    img(x+1, y+1, 2) = 0;
    img(x+1, y+1, 3) = 255;
end

%% outliers
outliers_img = outliers' * velo2rect * rect2disp;
outliers_img = outliers_img';
outliers_img = outliers_img ./ outliers_img(4, :);

outliers_visiable_idx = find(outliers_img(1, :) <= 1242-1 & ...
    outliers_img(1, :) >= 0 & ...
    outliers_img(2, :) <= 375-1 & ...
    outliers_img(2, :) >= 0 & ...
    outliers_img(3, :) <= 255 & ...
    outliers_img(3, :) >= 0);
outliers_img = outliers_img(:, outliers_visiable_idx);

outliers2d_img = round(outliers_img(1:2, :)) + 1;
for i = 1:size(outliers2d_img, 2)
    y = outliers2d_img(1, i);
    x = outliers2d_img(2, i);
    img(x, y, 1) = 255;
    img(x, y, 2) = 0;
    img(x, y, 3) = 0;
    
    img(x+1, y, 1) = 255;
    img(x+1, y, 2) = 0;
    img(x+1, y, 3) = 0;
    
    img(x, y+1, 1) = 255;
    img(x, y+1, 2) = 0;
    img(x, y+1, 3) = 0;
    
    img(x+1, y+1, 1) = 255;
    img(x+1, y+1, 2) = 0;
    img(x+1, y+1, 3) = 0;
end
Roadimage = img;
end

function [idInliers, idOutliers] = cut(X_tilde, B, th)
dists = abs(B'*normc(X_tilde));
idInliers = find(dists <= th);
idOutliers = find(dists > th);
end
