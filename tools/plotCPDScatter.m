function [fig1,fig2] = plotCPDScatter(popNeuronReg)
%plotCPDScatter Plot coefficients of partial determination during stimulus
%period (SV versus NSV or CV versus UCV)
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

x2 = r_cpd2(:,2);
y2 = r_cpd2(:,3);
sig_x2 = r_pValue2(:,2) <= threshold & r_pValue2(:,3) > threshold;
sig_y2 = r_pValue2(:,2) > threshold & r_pValue2(:,3) <= threshold;
sig_both2 = r_pValue2(:,2) <= threshold & r_pValue2(:,3) <= threshold;

% Plot figures
% Properties
max_lim = 0.45;
dotColor1 = [241,196,15;39,174,96;0,0,0]/255;
dotColor2 = [192,57,43;155,89,182;0,0,0]/255;
% Figure 1
fig1 = figure; hold on;
scatter(x1(sig_x1),y1(sig_x1),70,'filled','MarkerEdgeColor','k','MarkerFaceColor',dotColor1(1,:));
scatter(x1(sig_y1),y1(sig_y1),70,'filled','MarkerEdgeColor','k','MarkerFaceColor',dotColor1(2,:));
scatter(x1(sig_both1),y1(sig_both1),70,'filled','MarkerEdgeColor','k','MarkerFaceColor',dotColor1(3,:));
xlabel('CPD of salient value');
ylabel('CPD of non-salient value');
xlim([0 max_lim]); ylim([0 max_lim]);
line([0 max_lim],[0 max_lim],'lineStyle','--','color',[51 51 51]/255);
set(gca,'fontSize',13);
hold off
% Figure 2
fig2 = figure; hold on;
scatter(x2(sig_x2),y2(sig_x2),70,'filled','MarkerEdgeColor','k','MarkerFaceColor',dotColor2(1,:));
scatter(x2(sig_y2),y2(sig_y2),70,'filled','MarkerEdgeColor','k','MarkerFaceColor',dotColor2(2,:));
scatter(x2(sig_both2),y2(sig_both2),70,'filled','MarkerEdgeColor','k','MarkerFaceColor',dotColor2(3,:));
xlabel('CPD of cued value');
ylabel('CPD of uncued value');
xlim([0 max_lim]); ylim([0 max_lim]);
line([0 max_lim],[0 max_lim],'lineStyle','--','color',[51 51 51]/255);
set(gca,'fontSize',13);
hold off

end