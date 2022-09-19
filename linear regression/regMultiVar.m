function regResults = regMultiVar(popNeuron)

% Perform linear regression with regressors of attention cue location, SV
% and NSV on responses aligned to stimulus on and luminance change
fprintf('>>>> Performing linear regression on att_cue_loc, SV, and NSV ...\n');
temp_y = cellfun(@(x) x.norm_fr,popNeuron,'UniformOutput',false);
temp_x1 = cellfun(@(x) x.att_cue_loc,popNeuron,'UniformOutput',false);
temp_x2 = cellfun(@(x) x.SV,popNeuron,'UniformOutput',false);
temp_x3 = cellfun(@(x) x.NSV,popNeuron,'UniformOutput',false);
N = length(temp_y);
for i = 1:N
    fprintf(['>>>> Performing regression of neuron %d (',num2str(N),')... \n'],i);
    regResults.attCue_SV_NSV{i,1} ...
        = regSingleNeuron(temp_y{i},[temp_x1{i},temp_x2{i},temp_x3{i}],25);
end

% Perform linear regression with regressors of attention cue location, CV
% and UCV on responses aligned to stimulus on and luminance change
fprintf('>>>> Performing linear regression on att_cue_loc, CV, and UCV ...\n');
temp_y = cellfun(@(x) x.norm_fr,popNeuron,'UniformOutput',false);
temp_x1 = cellfun(@(x) x.att_cue_loc,popNeuron,'UniformOutput',false);
temp_x2 = cellfun(@(x) x.CV,popNeuron,'UniformOutput',false);
temp_x3 = cellfun(@(x) x.UCV,popNeuron,'UniformOutput',false);
N = length(temp_y);
for i = 1:N
    fprintf(['>>>> Performing regression of neuron %d (',num2str(N),')... \n'],i);
    regResults.attCue_CV_UCV{i,1} ...
        = regSingleNeuron(temp_y{i},[temp_x1{i},temp_x2{i},temp_x3{i}],25);
end

% Perform linear regression with regressors of attention cue location, SV
% and NSV on average responses during stimulus period
fprintf('>>>> Performing linear regression with average responses ...\n');
temp_y = cellfun(@(x) x.aveStim_norm_fr,popNeuron,'UniformOutput',false);
temp_x1 = cellfun(@(x) x.att_cue_loc,popNeuron,'UniformOutput',false);
temp_x2 = cellfun(@(x) x.SV,popNeuron,'UniformOutput',false);
temp_x3 = cellfun(@(x) x.NSV,popNeuron,'UniformOutput',false);
N = length(temp_y);
for i = 1:N
    regResults.aveStim_attCue_SV_NSV{i,1} ...
        = regSingleNeuron(temp_y{i},[temp_x1{i},temp_x2{i},temp_x3{i}]);
end

% Perform linear regression with regressors of attention cue location, CV
% and UCV on average responses during stimulus period
temp_y = cellfun(@(x) x.aveStim_norm_fr,popNeuron,'UniformOutput',false);
temp_x1 = cellfun(@(x) x.att_cue_loc,popNeuron,'UniformOutput',false);
temp_x2 = cellfun(@(x) x.CV,popNeuron,'UniformOutput',false);
temp_x3 = cellfun(@(x) x.UCV,popNeuron,'UniformOutput',false);
N = length(temp_y);
for i = 1:N
    regResults.aveStim_attCue_CV_UCV{i,1} ...
        = regSingleNeuron(temp_y{i},[temp_x1{i},temp_x2{i},temp_x3{i}]);
end

% Perform linear regression with regressors of attention cue location, SV
% and NSV on average responses after attention cue off
temp_y = cellfun(@(x) x.aveCueOff_norm_fr,popNeuron,'UniformOutput',false);
temp_x1 = cellfun(@(x) x.att_cue_loc,popNeuron,'UniformOutput',false);
temp_x2 = cellfun(@(x) x.SV,popNeuron,'UniformOutput',false);
temp_x3 = cellfun(@(x) x.NSV,popNeuron,'UniformOutput',false);
N = length(temp_y);
for i = 1:N
    regResults.aveCueOff_attCue_SV_NSV{i,1} ...
        = regSingleNeuron(temp_y{i},[temp_x1{i},temp_x2{i},temp_x3{i}]);
end

end

function regstats = regSingleNeuron(y, x, bin_width, varargin)

if nargin < 3 || isempty(bin_width), bin_width = 25; end

% Perform linear regression
if iscell(y)
    E = length(y);
    cpd = cell(1,E);
    coeff = cell(1,E);
    pValue = cell(1,E);
    for i = 1:E
        M = size(y{i},1);
        T = size(y{i},2)/bin_width;
        temp_y_bin = zeros(M,T);
        for j = 1:M
            temp_y_bin(j,:) = mean(reshape(y{i}(j,:),bin_width,[]));
        end        
        for k = 1:T
            [cpd{i}(:,k),coeff{i}(:,k),pValue{i}(:,k)] ...
                = regSingleDataPoint(temp_y_bin(:,k),x); 
        end
    end
    y_b = mean(y{1}(:,1:200),2);
    [cpd_b,~,~] = regSingleDataPoint(y_b,x);   
    regstats.cpd_b = cpd_b;
else
    [cpd,coeff,pValue] = regSingleDataPoint(y,x);    
end
regstats.cpd = cpd;
regstats.coeff = coeff;
regstats.pValue = pValue;
end