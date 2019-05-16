function [f1s] = f1score(error,label,ths)
    f1s = [];
    for th=ths
        assert(all(size(error) == size(label)), 'size of error and label should be equal');   
        assert(size(error, 1) == 1, 'error and label should be row array');
        inliers = find(label == 1);
        positive = find(error <= th);
        tp = intersect(inliers, positive);
        prec = length(tp) / length (positive);
        recall = length(tp) / length (inliers);
        f1 = 2*(prec*recall)/(prec+recall);
        if isnan(f1)
            f1 = 0;
        end
        f1s = [f1s; f1];
    end
end
