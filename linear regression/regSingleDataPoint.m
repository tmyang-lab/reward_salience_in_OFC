function [cpd,coeff,pValue] = regSingleDataPoint(y,x)

% Normalization
z = zscore(x,[],1);
% Perform linear regression
lm_full = fitlm(z,y);
col = size(z,2);
cpd = zeros(1,col);
for i = 1:col
    temp_z_red = z;
    temp_z_red(:,i) = [];
    temp_lm_red = fitlm(temp_z_red,y);
    cpd(i) = (temp_lm_red.SSE - lm_full.SSE)/temp_lm_red.SSE;
end
coeff = lm_full.Coefficients.Estimate(2:end)';
pValue = lm_full.Coefficients.pValue(2:end)';
end