function [bestB, best_inner, best_inliers, trialcount, timeAlgo] = RANSAC_PCA(X, Ns, threshold, budget)

[D, N] = size(X);

% parameters
p = 0.99;
maxTrials = 1000;

% initializations
T = maxTrials;
bestB = NaN;
best_inner = [];
trialcount = 0;
bestscore = 0;
t_start = tic;
while ((budget > toc(t_start)) && (T > trialcount))
    % Generate s random indices in the range 1..npts
    perm = randperm(N);
    ind = perm(1:Ns);
    
    % compute a subspace model
    B = fittingfn(X(:, ind), Ns);
    
    % Find the number of inliers to this model
    [inliers, inner, B] = distfn(B, X, threshold);
    ninliers = length(inliers);
    
    % update model if applicable
    if ninliers > bestscore
        bestscore = ninliers;
        best_inliers = inliers;
        bestB = B;
        best_inner = inner;
        
        % Update estimate of N, the number of trials to ensure we pick,
        % with probability p, a data set with no outliers.
        fracinliers = ninliers / N;
        pNoOutliers = 1 - fracinliers^Ns;
        pNoOutliers = max(eps, pNoOutliers);
        pNoOutliers = min(1-eps, pNoOutliers);
        T = min(log(1-p)/log(pNoOutliers), maxTrials);
    end
    
    trialcount = trialcount + 1;
end
timeAlgo = toc(t_start);

% disp(['RANSAC finished with ', num2str(trialcount), ' iterations.'])
