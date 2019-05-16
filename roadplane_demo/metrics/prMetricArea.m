function [AUCPR, AUCROC] = prMetricArea(errorArray, gtArray, vis)
if vis
    [prec, tpr, fpr, thresh] = prec_rec(-errorArray, gtArray, 'plotROC', 1, 'plotPR', 1);
else
    [prec, tpr, fpr, thresh] = prec_rec(-errorArray, gtArray);
end
AUCROC = trapz([0; fpr], [0; tpr]);
AUCPR = trapz([0; tpr], [1; prec]);
end
