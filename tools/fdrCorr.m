% This file is used to correct p value when performing multiple comparison.
% The correction method is fdr.
% input1: pValue is an 1 by m vector
% input2: alpha is a double
% output: hypothesis ia an 1 by m vector, 1 means alternative hypothesis is
% true, 0 means null hypothesis is true.
function hypothesis = fdrCorr(pValue,alpha)
% sort p value from the lowest to the highest, p_1 <= p_2 <= ... <= p_j
sort_pValue = sort(pValue);
% threshold is alpha*j/m, m is the length of pValue
m = length(pValue);
threshold = alpha*(1:1:m)/m;
% find the largest j for which p_j <= threshold
index = find(sort_pValue <= threshold,1,'last');
if isempty(index)
    hypothesis = zeros(1,m);
else
    % p_1, p_2, ..., p_j are significant if p_j <= threshold
    hypothesis = pValue <= sort_pValue(index);
end
end