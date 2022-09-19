function [fig1,fig2] = plotCPD2Areas(popNeuronReg1,popNeuronReg2)


% Extract variables
% Regression in area 1
N1 = length(popNeuronReg1);
cpd1 = cellfun(@(x) x.cpd,popNeuronReg1,'UniformOutput',false);
pValue1 = cellfun(@(x) x.pValue,popNeuronReg1,'UniformOutput',false);
% Regression in area 2
N2 = length(popNeuronReg2);
cpd2 = cellfun(@(x) x.cpd,popNeuronReg2,'UniformOutput',false);
pValue2 = cellfun(@(x) x.pValue,popNeuronReg2,'UniformOutput',false);

% Reshape variables
r_cpd1 = cell2mat(cpd1);
r_cpd2 = cell2mat(cpd2);
r_pValue1 = cell2mat(pValue1);
r_pValue2 = cell2mat(pValue2);

% Siginificant data points
threshold = 0.01;
x1 = r_cpd1(:,1);
y1 = r_cpd1(:,2);
x2 = r_cpd2(:,1);
y2 = r_cpd2(:,2);
sig1 = r_pValue1(:,1) <= threshold | r_pValue1(:,2) <= threshold;
sig2 = r_pValue2(:,1) <= threshold | r_pValue2(:,2) <= threshold;

% Calculate CPD difference
delta_cpd1 = x1 - y1;
delta_cpd2 = x2 - y2;
ave_sig_delta_cpd1 = mean(delta_cpd1(sig1));
ave_sig_delta_cpd2 = mean(delta_cpd2(sig2));

% Plot figures
% Properties
max_lim = 0.4;
dotColor = [255,95,31;0,255,255]/255;
% Figure 1
fig1 = figure; hold on;
line([0 max_lim],[0 max_lim],'lineStyle','--','color',[51 51 51]/255);
s2 = scatter(x2(sig2),y2(sig2),70,...
    'MarkerEdgeColor','k','MarkerFaceColor',dotColor(2,:),...
    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.8);
s1 = scatter(x1(sig1),y1(sig1),70,...
    'MarkerEdgeColor','k','MarkerFaceColor',dotColor(1,:),...
     'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.8);
legend([s1,s2],'OFC','LPFC');
legend('boxoff');
xlabel('CPD of attention cue location');
ylabel('CPD of salient value');
xlim([0 max_lim]); ylim([0 max_lim]);
hold off

% Figure 2
fig2 = figure('Position',[680,558,560,230]); hold on;
edges = -max_lim:max_lim/25:max_lim;
edges_bin = -max_lim+max_lim/25/2:max_lim/25:max_lim-max_lim/25/2;
line([0 0],[0 1],'lineStyle','--','color',[51 51 51]/255);
N_delta_cpd2 = histcounts(delta_cpd2(sig2),edges);
bar(edges_bin,N_delta_cpd2/N2,'EdgeColor','k','FaceColor',dotColor(2,:),...
    'EdgeAlpha',1,'FaceAlpha',0.5);
N_delta_cpd1 = histcounts(delta_cpd1(sig1),edges);
bar(edges_bin,N_delta_cpd1/N1,'EdgeColor','k','FaceColor',dotColor(1,:),...
    'EdgeAlpha',1,'FaceAlpha',0.5);
line([ave_sig_delta_cpd2 ave_sig_delta_cpd2],[0 0.15],...
    'lineStyle','--','color',dotColor(2,:),'lineWidth',2);
line([ave_sig_delta_cpd1 ave_sig_delta_cpd1],[0 0.15],...
    'lineStyle','--','color',dotColor(1,:),'lineWidth',2);
xlabel('Delta CPD');
ylabel('Proportion of units');
xlim([-max_lim max_lim]); ylim([0 0.06]);
hold off

end