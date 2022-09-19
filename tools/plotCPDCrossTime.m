function [fig1,fig2] = plotCPDCrossTime(popNeuronReg)
%plotCPDCrossTime Plot coefficients of partial determination at each time
%bin
% Input: 
%    popNeuronReg: is a structure, containing regression results of each
%    neruon
% Output:
%    fig1: CPD of regression with regressors of attention cue location, SV,
%    and NSV
%    fig2: CPD of regression with regressors of attention cue location, CV
%    and UCV
    
% Extract variables
% Regression with regressors of attention cue location, SV and NSV
cpd_b1 = cellfun(@(x) x.cpd_b,popNeuronReg.attCue_SV_NSV,'UniformOutput',false);
cpd1 = cellfun(@(x) x.cpd,popNeuronReg.attCue_SV_NSV,'UniformOutput',false);
pValue1 = cellfun(@(x) x.pValue,popNeuronReg.attCue_SV_NSV,'UniformOutput',false);
% Regression with regressors of attention cue location, CV and UCV
cpd_b2 = cellfun(@(x) x.cpd_b,popNeuronReg.attCue_CV_UCV,'UniformOutput',false);
cpd2 = cellfun(@(x) x.cpd,popNeuronReg.attCue_CV_UCV,'UniformOutput',false);
pValue2 = cellfun(@(x) x.pValue,popNeuronReg.attCue_CV_UCV,'UniformOutput',false);

% Reshape variables
r_cpd_b1 = cell2mat(cpd_b1);
r_cpd1 = cell(1,2);
r_pValue1 = cell(1,2);
r_cpd_b2 = cell2mat(cpd_b2);
r_cpd2 = cell(1,2);
r_pValue2 = cell(1,2);
N = length(cpd1);
for i = 1:N
    for j = 1:2
        r_cpd1{j}(i,:,:) = cpd1{i}{j};
        r_pValue1{j}(i,:,:) = pValue1{i}{j};
        r_cpd2{j}(i,:,:) = cpd2{i}{j};
        r_pValue2{j}(i,:,:) = pValue2{i}{j};
    end
end

% Significance test
hypothesis1 = cell(1,2);
hypothesis2 = cell(1,2);
for i = 1:2
    T = size(r_cpd1{i},3);
    temp_var_cpd1 = zeros(N,3,T);
    temp_p1 = zeros(4,T);
    temp_var_cpd2 = zeros(N,3,T);
    temp_p2 = zeros(4,T);
    for j = 1:T
        for k = 1:3
            temp_var_cpd1(:,k,j) = reshape(r_cpd1{i}(:,k,j),[N,1]);
            [~,temp_p1(k,j)] = ttest(temp_var_cpd1(:,k,j),mean(r_cpd_b1,2));
            temp_var_cpd2(:,k,j) = reshape(r_cpd2{i}(:,k,j),[N,1]);
            [~,temp_p2(k,j)] = ttest(temp_var_cpd2(:,k,j),mean(r_cpd_b2,2));
        end
        [~,temp_p1(4,j)] = ttest(temp_var_cpd1(:,2,j),temp_var_cpd1(:,3,j));
        [~,temp_p2(4,j)] = ttest(temp_var_cpd2(:,2,j),temp_var_cpd2(:,3,j));
    end
    hypothesis1{i} = reshape(fdrCorr(temp_p1(:)',0.005),[4,T]);
    hypothesis2{i} = reshape(fdrCorr(temp_p2(:)',0.005),[4,T]);
end

% Plot figures
% Properties
time_label{1} = 'Time from stim on (ms)';
time_label{2} = 'Time from lum change (ms)';
tw{1} = [-400 600];
tw{2} = [-400 400];
bin_width = 25;
shadeColor1 = [41,128,185;241,196,15;39,174,96;0,0,0]/255;
shadeColor2 = [41,128,185;192,57,43;155,89,182;0,0,0]/255;
% Figure 1
fig1 = figure('Position',[500 500 800 400]);
for i = 1:2  
    subplot(1,2,i); hold on;
    T = size(r_cpd1{i},3);
    X = tw{i}(1)+bin_width/2:bin_width:tw{i}(end)-bin_width/2;
    y = [ones(1,T)*0.041;ones(1,T)*0.042;ones(1,T)*0.043;ones(1,T)*0.044]; 
    for j = 1:4
        scatter(X(hypothesis1{i}(j,:)==1),y(j,hypothesis1{i}(j,:)==1),...
        'filled','s','MarkerEdgeColor',shadeColor1(j,:),'MarkerFaceColor',shadeColor1(j,:));
    end
    for j = 1:3
        temp_ave_cpd = reshape(mean(r_cpd1{i}(:,j,:),1),[1,T]);
        temp_sem_cpd = reshape(std(r_cpd1{i}(:,j,:),1)/sqrt(N),[1,T]);
        e(j) = errorbar(X,temp_ave_cpd,temp_sem_cpd,'color',shadeColor1(j,:),'lineWidth',2);
    end  
    if i == 2
        legend([e(1),e(2),e(3)],'Attention cue loc','Salient value','Non-salient value','Location','West');
        legend('boxoff');
    end
    xlim([tw{i}(1) tw{i}(end)]); 
    ylim([0 0.045]);
    xlabel(time_label{i}); ylabel('CPD');
    set(gca,'fontSize',13);
    hold off;
end

% Figure 2
fig2 = figure('Position',[500 500 800 400]);
for i = 1:2  
    subplot(1,2,i); hold on;
    T = size(r_cpd2{i},3);
    X = tw{i}(1)+bin_width/2:bin_width:tw{i}(end)-bin_width/2;
    y = [ones(1,T)*0.041;ones(1,T)*0.042;ones(1,T)*0.043;ones(1,T)*0.044];  
    for j = 1:4
        scatter(X(hypothesis2{i}(j,:)==1),y(j,hypothesis2{i}(j,:)==1),...
        'filled','s','MarkerEdgeColor',shadeColor2(j,:),'MarkerFaceColor',shadeColor2(j,:));
    end
    for j = 1:3
        temp_ave_cpd = reshape(mean(r_cpd2{i}(:,j,:),1),[1,T]);
        temp_sem_cpd = reshape(std(r_cpd2{i}(:,j,:),1)/sqrt(N),[1,T]);
        e(j) = errorbar(X,temp_ave_cpd,temp_sem_cpd,'color',shadeColor2(j,:),'lineWidth',2);
    end  
    if i == 2
        legend([e(1),e(2),e(3)],'Attention cue loc','Cued value','Uncued value','Location','West');
        legend('boxoff');
    end
    xlim([tw{i}(1) tw{i}(end)]); 
    ylim([0 0.045]);
    xlabel(time_label{i}); ylabel('CPD');
    set(gca,'fontSize',13);
    hold off;
end

end