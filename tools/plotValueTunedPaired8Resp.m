function fig = plotValueTunedPaired8Resp(popNeuron)
%plotValueTunedPaired8Resp Plot responses of population OFC neurons under
%values paired with 8-drops-of-juice stimulus
% Inputs: 
%    popNeuron: is a neuron num by 1 cell array, each cell contains a
%    structure, containing responses and task-related variables of the
%    neuron
% Output: 
%    fig: firing rates sorted by same value, UCVs paired with cued
%    8-drops-of-juice stimulus, CVs paired with uncued-8-drops-of juice
%    stimulus.

% Extract variables
aveStim_fr = cellfun(@(x) x.aveStim_norm_fr,popNeuron,'UniformOutput',false);
SaV = cellfun(@(x) x.SaV,popNeuron,'UniformOutput',false);
CV = cellfun(@(x) x.CV,popNeuron,'UniformOutput',false);
UCV = cellfun(@(x) x.UCV,popNeuron,'UniformOutput',false);
% Sort trials
value_cond = [0 1 2 4 8];
aveStim_fr_same = cell(1,5);
aveStim_fr_cued8 = cell(1,5);
aveStim_fr_uncued8 = cell(1,5);
for i = 1:length(aveStim_fr)
    for j = 1:5
        aveStim_fr_same{j}(i,1) = mean(aveStim_fr{i}(SaV{i} == value_cond(j))); 
        aveStim_fr_cued8{j}(i,1) = mean(aveStim_fr{i}(CV{i} == 8 & UCV{i} == value_cond(j)));
        aveStim_fr_uncued8{j}(i,1) = mean(aveStim_fr{i}(UCV{i} == 8 & CV{i} == value_cond(j)));
    end
end

% Plot figure
% Properties
lineColor = [0,0,0;192,57,43;155,89,182]/255;
% Figure 
fig = figure; hold on;
temp_fr1 = aveStim_fr_same;
temp_fr2 = aveStim_fr_cued8;
temp_fr3 = aveStim_fr_uncued8;
ave_temp_fr1 = zeros(1,5);
sem_temp_fr1 = zeros(1,5);
ave_temp_fr2 = zeros(1,5);
sem_temp_fr2 = zeros(1,5);
ave_temp_fr3 = zeros(1,5);
sem_temp_fr3 = zeros(1,5);
for i = 1:5
    ave_temp_fr1(1,i) = nanmean(temp_fr1{i},1);
    sem_temp_fr1(1,i) = nanstd(temp_fr1{i},0,1)/sqrt(length(temp_fr1{i}));
    ave_temp_fr2(1,i) = nanmean(temp_fr2{i},1);
    sem_temp_fr2(1,i) = nanstd(temp_fr2{i},0,1)/sqrt(length(temp_fr2{i}));
    ave_temp_fr3(1,i) = nanmean(temp_fr3{i},1);
    sem_temp_fr3(1,i) = nanstd(temp_fr3{i},0,1)/sqrt(length(temp_fr3{i}));
end
e2 = errorbar(value_cond,ave_temp_fr2,sem_temp_fr2,'Color',lineColor(2,:),'LineWidth',2);
e3 = errorbar(value_cond,ave_temp_fr3,sem_temp_fr3,'Color',lineColor(3,:),'LineWidth',2);
e1 = errorbar(value_cond,ave_temp_fr1,sem_temp_fr1,'Color',lineColor(1,:),'LineWidth',2);
legend([e1,e2,e3],{'Same value','Cued value = 8','Uncued value = 8'},'Location','North');
legend('boxoff');
xlabel('Value size'); 
ylabel('Normalized firing rates (spikes/s)');
xlim([-1,8]); 
ylim([-0.6,0.6]);
set(gca,'FontSize',14);
hold off;

end