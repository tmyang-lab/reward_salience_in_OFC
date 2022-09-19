function fig = plotExampleDLPFCResp(expNeuron)
%plotExampleDLPFCResp Plot responses of the example DLPFC neuron
% Inputs: 
% Inputs: 
%    expNeuron: is a structure, containing responses and task-related 
%    variables of the neuron
% Output: 
%    fig: firing rates sorted by attention cue location

% Extract variables
y = expNeuron.smooth_fr;
x = expNeuron.att_cue_loc;
% Sort trials 
y_sort{1}{1} = y{1}(x == 1,:);
y_sort{1}{2} = y{1}(x == -1,:);
y_sort{2}{1} = y{2}(x == 1,:);
y_sort{2}{2} = y{2}(x == -1,:);

% Plot figure
% Properties
timeLabel{1} = 'Time from stim on (ms)';
timeLabel{2} = 'Time from lum change (ms)'; 
tw{1} = [-400 600];
tw{2} = [-400 400];
line_color = [52 152 219;231 76 60]/255;
fig = figure('Position',[500,500,800,400]);
for i = 1:2    
    subplot(1,2,i);hold on;
    X = tw{i}(1)+0.5:1:tw{i}(end)-0.5;  
    for j = 1:2
        temp_ave = nanmean(y_sort{i}{j},1);
        temp_sem = nanstd(y_sort{i}{j},0,1)/sqrt(size(y_sort{i}{j},1));       
        patch([X fliplr(X)],[temp_ave+temp_sem ...
            fliplr(temp_ave-temp_sem)],...
            [160 160 160]/255,...
            'facealpha',0.5,...
            'EdgeColor',[160 160 160]/255,...
            'edgealpha',0.5);
    end
    for j = 1:2
        temp_ave = nanmean(y_sort{i}{j},1);
        p(j) = plot(X,temp_ave,'color',line_color(j,:),'lineWidth',2); 
    end 
    if i == 1
        legend([p(2),p(1)],...
            'Left attention cue','Right attention cue',...
            'Location','northeast');
        legend('boxoff');
    end
    xlim([tw{i}(1) tw{i}(end)]);
    ylim([0 35]);
    xlabel(timeLabel{i}); 
    ylabel('Firing rates (spikes/s)');
%     set(gca,'FontSize',12);
    hold off;
end
end