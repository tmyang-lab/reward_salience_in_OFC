function [fig1,fig2,fig3,fig4] = plotValueTunedResp(popNeuron)
%plotValueTunedResp Plot responses of the population OFC neuron under
%different values
% Inputs: 
%    popNeuron: is a neuron num by 1 cell array, each cell contains a
%    structure, containing responses and task-related variables of the
%    neuron
% Output: 
%    -(1) fig1: firing rates sorted by same value
%    -(2) fig2: firing rates sorted by salient value
%    -(3) fig3: compare firing rates sorted by same value and salient value
%    -(4) fig4: compare firing rates sorted by cued value and un-cued value

% Extract variables
fr = cellfun(@(x) x.norm_fr,popNeuron,'UniformOutput',false);
aveStim_fr = cellfun(@(x) x.aveStim_norm_fr,popNeuron,'UniformOutput',false);
SaV = cellfun(@(x) x.SaV,popNeuron,'UniformOutput',false);
SV = cellfun(@(x) x.SV,popNeuron,'UniformOutput',false);
CV = cellfun(@(x) x.CV,popNeuron,'UniformOutput',false);
UCV = cellfun(@(x) x.UCV,popNeuron,'UniformOutput',false);
% Sort trials
value_cond = [0 1 2 4 8];
fr_same = cell(1,2);
fr_salient = cell(1,2);
fr_cued = cell(1,2);
fr_uncued = cell(1,2);
for i = 1:length(fr)
    for j = 1:2
        for k = 1:5
            fr_same{j}{k}(i,:) = mean(fr{i}{j}(SaV{i} == value_cond(k),:),1); 
            fr_salient{j}{k}(i,:) = mean(fr{i}{j}(SV{i} == value_cond(k),:),1);
            fr_cued{j}{k}(i,:) = mean(fr{i}{j}(CV{i} == value_cond(k),:),1);
            fr_uncued{j}{k}(i,:) = mean(fr{i}{j}(UCV{i} == value_cond(k),:),1);        
        end
    end
end
aveStim_fr_same = cell(1,5);
aveStim_fr_salient = cell(1,5);
aveStim_fr_cued = cell(1,5);
aveStim_fr_uncued = cell(1,5);
for i = 1:length(aveStim_fr)
    for j = 1:5
        aveStim_fr_same{j}(i,1) = mean(aveStim_fr{i}(SaV{i} == value_cond(j))); 
        aveStim_fr_salient{j}(i,1) = mean(aveStim_fr{i}(SV{i} == value_cond(j)));
        aveStim_fr_cued{j}(i,1) = mean(aveStim_fr{i}(CV{i} == value_cond(j)));
        aveStim_fr_uncued{j}(i,1) = mean(aveStim_fr{i}(UCV{i} == value_cond(j)));
    end
end

% Plot figure
% Properties
timeLabel{1} = 'Time from stim on (ms)';
timeLabel{2} = 'Time from lum change (ms)'; 
tw{1} = [-400 600];
tw{2} = [-400 400];
shadeColor1 = [224 224 224;192 192 192;128 128 128;64 64 64;0 0 0]/255;
shadeColor2 = [249 231 159;247 220 111;241 196 15;183 149 11;125 102 8]/255;
lineColor1 = [0 0 0;241 196 15]/255;
lineColor2 = [192 57 43;155 89 182]/255;
bin_width = 1; 

% Figure1
fig1 = figure('Position',[500,500,800,400]);
temp_fr = fr_same;
for i = 1:2    
    subplot(1,2,i);hold on;
    X = tw{i}(1)+bin_width/2:bin_width:tw{i}(end)-bin_width/2;  
    for j = 1:5        
        ave_fr = nanmean(temp_fr{i}{j},1);
        sem_fr = nanstd(temp_fr{i}{j},0,1)/sqrt(size(temp_fr{i}{j},1));
        patch([X fliplr(X)],[ave_fr+sem_fr fliplr(ave_fr-sem_fr)],...
            [160 160 160]/255,'facealpha',0.5,'EdgeColor',[160 160 160]/255,'edgealpha',0.5);
    end
    for j = 1:5
        ave_fr = nanmean(temp_fr{i}{j},1);
        p(j) = plot(X,ave_fr,'color',shadeColor1(j,:),'lineWidth',2); 
    end 
    if i == 1
        legend([p(1),p(2),p(3),p(4),p(5)],...
            {'~0 drop','1 drop','2 drops','4 drops','8 drops'},'Location','NorthWest');
        legend('boxoff');
    end
    xlim([tw{i}(1) tw{i}(end)]);
    ylim([-0.5 0.6]);
    xlabel(timeLabel{i},'FontSize',12); ylabel('Normalized firing rates (spikes/s)','FontSize',12);
    set(gca,'FontSize',12);
    hold off;
end

% Figure 2
fig2 = figure('Position',[500,500,800,400]);
temp_fr = fr_salient;
for i = 1:2    
    subplot(1,2,i);hold on;
    X = tw{i}(1)+bin_width/2:bin_width:tw{i}(end)-bin_width/2;  
    for j = 1:5        
        ave_fr = nanmean(temp_fr{i}{j},1);
        sem_fr = nanstd(temp_fr{i}{j},0,1)/sqrt(size(temp_fr{i}{j},1));
        patch([X fliplr(X)],[ave_fr+sem_fr fliplr(ave_fr-sem_fr)],...
            [160 160 160]/255,'facealpha',0.5,'EdgeColor',[160 160 160]/255,'edgealpha',0.5);
    end
    for j = 1:5
        ave_fr = nanmean(temp_fr{i}{j},1);
        p(j) = plot(X,ave_fr,'color',shadeColor2(j,:),'lineWidth',2); 
    end 
    if i == 1
        legend([p(1),p(2),p(3),p(4),p(5)],...
            {'~0 drop','1 drop','2 drops','4 drops','8 drops'},'Location','NorthWest');
        legend('boxoff');
    end
    xlim([tw{i}(1) tw{i}(end)]);
    ylim([-0.5 0.6]);
    xlabel(timeLabel{i},'FontSize',12); ylabel('Normalized firing rates (spikes/s)','FontSize',12);
    set(gca,'FontSize',12);
    hold off;
end

% Figure 3
fig3 = figure; hold on;
temp_fr1 = aveStim_fr_same;
temp_fr2 = aveStim_fr_salient;
ave_temp_fr1 = zeros(1,5);
sem_temp_fr1 = zeros(1,5);
ave_temp_fr2 = zeros(1,5);
sem_temp_fr2 = zeros(1,5);
for i = 1:5
    ave_temp_fr1(1,i) = nanmean(temp_fr1{i},1);
    sem_temp_fr1(1,i) = nanstd(temp_fr1{i},0,1)/sqrt(length(temp_fr1{i}));
    ave_temp_fr2(1,i) = nanmean(temp_fr2{i},1);
    sem_temp_fr2(1,i) = nanstd(temp_fr2{i},0,1)/sqrt(length(temp_fr2{i}));
end
e2 = errorbar(value_cond,ave_temp_fr2,sem_temp_fr2,'Color',lineColor1(2,:),'LineWidth',2);
e1 = errorbar(value_cond,ave_temp_fr1,sem_temp_fr1,'Color',lineColor1(1,:),'LineWidth',2);
legend([e1,e2],{'Same value','Salient value'});
legend('boxoff');
xlabel('Value size'); ylabel('Normalized firing rates (spikes/s)');
xlim([-1,8]); 
ylim([-0.5 0.6]);
set(gca,'FontSize',14);
hold off;

% Figure 4
fig4 = figure; hold on;
temp_fr1 = aveStim_fr_cued;
temp_fr2 = aveStim_fr_uncued;
ave_temp_fr1 = zeros(1,5);
sem_temp_fr1 = zeros(1,5);
ave_temp_fr2 = zeros(1,5);
sem_temp_fr2 = zeros(1,5);
for i = 1:5
    ave_temp_fr1(1,i) = nanmean(temp_fr1{i},1);
    sem_temp_fr1(1,i) = nanstd(temp_fr1{i},0,1)/sqrt(length(temp_fr1{i}));
    ave_temp_fr2(1,i) = nanmean(temp_fr2{i},1);
    sem_temp_fr2(1,i) = nanstd(temp_fr2{i},0,1)/sqrt(length(temp_fr2{i}));
end
e1 = errorbar(value_cond,ave_temp_fr1,sem_temp_fr1,'Color',lineColor2(1,:),'LineWidth',2);
e2 = errorbar(value_cond,ave_temp_fr2,sem_temp_fr2,'Color',lineColor2(2,:),'LineWidth',2);
legend([e1,e2],'Cued value','Uncued value');
legend('boxoff');
xlabel('Value size'); ylabel('Normalized firing rates(spikes/s)');
xlim([-1,8]); 
ylim([-0.5 0.6]);
set(gca,'FontSize',14);
hold off;

end