function [fig1,fig2] = plotCPDHistogrm(popNeuronReg)
%plotCPDHistogrm Plot distribution of coefficients of partial determination 
%during stimulus period
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
cpd1 = cellfun(@(x) x.cpd,popNeuronReg.aveStim_attCue_SV_NSV,'UniformOutput',false);
pValue1 = cellfun(@(x) x.pValue,popNeuronReg.aveStim_attCue_SV_NSV,'UniformOutput',false);
% Regression with regressors of attention cue location, CV and UCV
cpd2 = cellfun(@(x) x.cpd,popNeuronReg.aveStim_attCue_CV_UCV,'UniformOutput',false);
pValue2 = cellfun(@(x) x.pValue,popNeuronReg.aveStim_attCue_CV_UCV,'UniformOutput',false);

% Reshape variables
r_cpd1 = cell2mat(cpd1);
r_cpd2 = cell2mat(cpd2);
r_pValue1 = cell2mat(pValue1);
r_pValue2 = cell2mat(pValue2);

% Siginificant data points
threshold = 0.05;
x1 = r_cpd1(:,2);
y1 = r_cpd1(:,3);
sig_x1 = r_pValue1(:,2) <= threshold & r_pValue1(:,3) > threshold;
sig_y1 = r_pValue1(:,2) > threshold & r_pValue1(:,3) <= threshold;
sig_both1 = r_pValue1(:,2) <= threshold & r_pValue1(:,3) <= threshold;
ave_x1 = nanmean(x1(sig_x1|sig_y1|sig_both1));
ave_y1 = nanmean(y1(sig_x1|sig_y1|sig_both1));

x2 = r_cpd2(:,2);
y2 = r_cpd2(:,3);
sig_x2 = r_pValue2(:,2) <= threshold & r_pValue2(:,3) > threshold;
sig_y2 = r_pValue2(:,2) > threshold & r_pValue2(:,3) <= threshold;
sig_both2 = r_pValue2(:,2) <= threshold & r_pValue2(:,3) <= threshold;
ave_x2 = nanmean(x2(sig_x2|sig_y2|sig_both2));
ave_y2 = nanmean(y2(sig_x2|sig_y2|sig_both2));

% Plot figures
% Properties
max_lim = 0.4;
barColor11 = [255,255,255;125,102,8]/255;
barColor12 = [255,255,255;20,90,50]/255;
barColor21 = [255,255,255;100,30,22]/255;
barColor22 = [255,255,255;81,46,95]/255;
% Figure 1
fig1 = figure; 
subplot(2,1,1); hold on;
edges = 0:max_lim/20:max_lim;
h2 = histogram(x1(sig_x1|sig_both1|sig_y1),edges,'EdgeColor','k','FaceColor',barColor11(1,:),'FaceAlpha',0.8);
h1 = histogram(x1(sig_x1|sig_both1),edges,'EdgeColor','k','FaceColor',barColor11(2,:),'FaceAlpha',0.8);    
line([ave_x1 ave_x1],[0 100],'lineStyle','--','color','k');
legend([h1,h2],'Sig salient value','Non-sig salient value');
legend('boxoff');
xlabel('CPD of salient value');
ylabel('Number of units');
xlim([-max_lim/10 max_lim]); 
% ylim([0 50]);
hold off;
subplot(2,1,2); hold on;
edges = 0:max_lim/20:max_lim;
h2 = histogram(y1(sig_y1|sig_both1|sig_x1),edges,'EdgeColor','k','FaceColor',barColor12(1,:),'FaceAlpha',0.8);
h1 = histogram(y1(sig_y1|sig_both1),edges,'EdgeColor','k','FaceColor',barColor12(2,:),'FaceAlpha',0.8);    
line([ave_y1 ave_y1],[0 100],'lineStyle','--','color','k');
legend([h1,h2],'Sig non-salient value','Non-sig non-salient value');
legend('boxoff');
xlabel('CPD of non-salient value');
ylabel('Number of units');
xlim([-max_lim/10 max_lim]); 
% ylim([0 100]);
hold off;

% Figure 2
fig2 = figure; 
subplot(2,1,1); hold on;
edges = 0:max_lim/20:max_lim;
h2 = histogram(x2(sig_x2|sig_both2|sig_y2),edges,'EdgeColor','k','FaceColor',barColor21(1,:),'FaceAlpha',0.8);
h1 = histogram(x2(sig_x2|sig_both2),edges,'EdgeColor','k','FaceColor',barColor21(2,:),'FaceAlpha',0.8);    
line([ave_x2 ave_x2],[0 100],'lineStyle','--','color','k');
legend([h1,h2],'Sig cued value','Non-sig cued value');
legend('boxoff');
xlabel('CPD of salient value');
ylabel('Number of units');
xlim([-max_lim/10 max_lim]); 
% ylim([0 50]);
hold off;
subplot(2,1,2); hold on;
edges = 0:max_lim/20:max_lim;
h2 = histogram(y2(sig_y2|sig_both2|sig_x2),edges,'EdgeColor','k','FaceColor',barColor22(1,:),'FaceAlpha',0.8);
h1 = histogram(y2(sig_y2|sig_both2),edges,'EdgeColor','k','FaceColor',barColor22(2,:),'FaceAlpha',0.8);    
line([ave_y2 ave_y2],[0 100],'lineStyle','--','color','k');
legend([h1,h2],'Sig uncued value','Non-sig uncued value');
legend('boxoff');
xlabel('CPD of non-salient value');
ylabel('Number of units');
xlim([-max_lim/10 max_lim]); 
% ylim([0 100]);
hold off;

end