clear all;
addAllPaths;

%% Setup
% specify the frame that we would like to work on
date = 29;
seq = 71;
frame = 328;

% switch for 3d point cloud visualization
vis_pc = 0;

% number of trials to average RANSAC
% note: nExp=1 is for getting a quick idea. for rigorous benchmark, one
% should average both RANSAC and other methods for some large number of trials.
nExp = 1;


%% Load and Prepare Data
% load the dataset
dataPath = './data/';
f = load([dataPath, sprintf('09_%02d_%03d/%02d_%03d_%04d_data.mat', date, seq, date, seq, frame)]);

% X, O: inliers, outliers in R^4 using homogenized coordinates
X = normc(f.inliers);
O = normc(f.outliers);
X_tilde = [X, O];

% # inliers, # outliers and inlier ratio
n_points = length(X_tilde);
n_inliers = size(X, 2);
n_outliers = size(O, 2);
inlier_ratio = n_inliers / n_points;

% produce ground truth normal in R^4
[V, d] = eig(X*X');
[d, ind] = sort(diag(d), 'descend');
V = V(:, ind);
B_gt = V(:, end);

% ground truth threshold
th_gt = max(abs(B_gt'*X));

% point clouds in R^3
points3D = X_tilde(1:3, :) ./ X_tilde(4, :);

% ground truth labeling
L = [ones(1, size(X, 2)), zeros(1, size(O, 2))];

%% run DPCP-PSGM
beta = 1 / 2;
[t_PSGM, B_PSGM, ~, iter_PSGM] = DPCP_PSGM_optim_J(X_tilde, 1, 1e-9, beta, 50, 1e-6);
[resultPSGM] = metric(B_PSGM, X_tilde, L, B_gt, th_gt, points3D, vis_pc);
resultPSGM.t = t_PSGM;
resultPSGM.iter = iter_PSGM

%% run DPCP-d
[B_d, t_d, iter_d] = DPCA_d(X_tilde, 100, 1e-6, 1, 2.8/sqrt(n_points), 1, inf);
[resultD] = metric(B_d, X_tilde, L, B_gt, th_gt, points3D, vis_pc);
resultD.t = t_d;
resultD.iter = iter_d

%% run REAPER-IRLS
[t_REAPER, ~, B_REAPER, iter_REAPER] = REAPER_IRLS_optim(X_tilde, 1, 1e-6, 100, 1e-9, t_PSGM);
[resultREAPER] = metric(B_REAPER, X_tilde, L, B_gt, th_gt, points3D, vis_pc);
resultREAPER.t = t_REAPER;
resultREAPER.iter = iter_REAPER

%% run DPCP-IRLS
[B_IRLS, t_IRLS, iter_IRLS] = DPCP_IRLS(X_tilde, 1, 1e-9, 100, 1e-6, t_PSGM);
[resultIRLS] = metric(B_IRLS, X_tilde, L, B_gt, th_gt, points3D, vis_pc);
resultIRLS.t = t_IRLS;
resultIRLS.iter = iter_IRLS

%% run L21-RPCA
[t_L21, ~, B_L21, iter_L21] = RPCA_L21_OutlierDetection(X_tilde, 1, 4*3/(7 * sqrt(n_outliers)), t_PSGM);
[resultL21] = metric(B_L21, X_tilde, L, B_gt, th_gt, points3D, vis_pc);
resultL21.t = t_L21;
resultL21.iter = iter_L21

%% run svd
[V, d] = eig(X_tilde*X_tilde');
[d, ind] = sort(diag(d), 'descend');
V = V(:, ind);
B_svd = V(:, end);
[resultSVD] = metric(B_svd, X_tilde, L, B_gt, th_gt, points3D, vis_pc);
resultSVD.t = 0;
resultSVD.iter = 1;

%% run RANSAC
% note: this takes a moderate time, since in order to plot a ROC curve, one
% has to vary the internal thresholding parameters, run RANSAC and repeat.
[resultRANSAC] = runRansac(X_tilde, th_gt, B_gt, L, t_PSGM, nExp, points3D, vis_pc)
[resultRANSAC_10] = runRansac(X_tilde, th_gt, B_gt, L, 10*t_PSGM, nExp, points3D, vis_pc)
[resultRANSAC_100] = runRansac(X_tilde, th_gt, B_gt, L, 100*t_PSGM, nExp, points3D, vis_pc)

%% visualization
resultsAll = struct();
resultsAll.PSGM = resultPSGM;
resultsAll.D = resultD;
resultsAll.REAPER = resultREAPER;
resultsAll.IRLS = resultIRLS;
resultsAll.L21 = resultL21;
resultsAll.SVD = resultSVD;
resultsAll.RANSAC = resultRANSAC;
resultsAll.RANSAC_10 = resultRANSAC_10;
resultsAll.RANSAC_100 = resultRANSAC_100;
resultsAll.n_inliers = n_inliers;
resultsAll.n_points = n_points;
resultsAll.n_outliers = n_outliers;
resultsAll.inlier_ratio = inlier_ratio;

visualizerIMG(date, seq, frame, dataPath, resultsAll);