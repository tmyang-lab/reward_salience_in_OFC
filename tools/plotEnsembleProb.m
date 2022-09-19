function [fig1,fig2] = plotEnsembleProb(ensembleNeuron)
%plotEnsembleProb Plot shaded area of posterior probability 

% Input: 
%    ensembleNeuron: is a structure, containing LDA results of neuronal
%    ensemble
% Output:
%    fig1: posterior probability from the decoder trained and tested with  
%    the responses in the same time bin
%    fig2: posterior probability from the decoder trained with the 
%    responses before stimulus onset

% Extract variables
P = ensembleNeuron.P;
shuf_P = ensembleNeuron.shuf_P;
P_f = ensembleNeuron.P_f;
shuf_P_f = ensembleNeuron.shuf_P_f;
% Properties
% Bin width
bin_width = 25;
% Step 
bin_step = 10;
tw{1} = [-400 600];
tw{2} = [-400 400];
timeLabel{1} = 'Time from stim on (ms)';
timeLabel{2} = 'Time from lum change (ms)';
% There are four matrix, each row is the training time, col is the testing
% time
% X is the time of each epoch align to imgage on and luminance change
% X_num is the length of each time epoch
epoch_num = 2;
X = cell(1,epoch_num);
X_num = zeros(1,epoch_num);
for i = 1:epoch_num
    X_ini = 1+tw{i}(1):bin_step:tw{i}(end)-bin_width+1;
    X_end = tw{i}(1)+bin_width:bin_step:tw{i}(end);
    X{i} = (X_ini + X_end)/2;
    X_num(i) = length(X{i});
end
X_num_extend = [0 cumsum(X_num)];
% X_index is the index of POSTERIOR corresponding to each time epoch
X_index = cell(1,epoch_num);
for i = 1:epoch_num
    X_index{i} = X_num_extend(i)+1:1:X_num_extend(i+1);
end
POSTERIOR = cell(1,epoch_num);
shuf_POSTERIOR = cell(1,epoch_num);
for i = 1:epoch_num
        POSTERIOR{i} = P(:,X_index{i});
        shuf_POSTERIOR{i} = shuf_P(:,X_index{i});
end
POSTERIOR_f = cell(1,epoch_num);
shuf_POSTERIOR_f = cell(1,epoch_num);
for i = 1:epoch_num
        POSTERIOR_f{i} = P_f(:,X_index{i});
        shuf_POSTERIOR_f{i} = shuf_P_f(:,X_index{i});
end

% Plot figure
fig1 = figure('Position',[500 500 800 400]);
for i = 1:epoch_num    
    subplot(1,epoch_num,i);hold on;
    M= size(P,1);
    ave_POSTERIOR = mean(POSTERIOR{i},1);
    sem_POSTERIOR = std(POSTERIOR{i},0,1)/sqrt(M);
    ave_shuf_POSTERIOR = mean(shuf_POSTERIOR{i},1);
    sem_shuf_POSTERIOR = std(shuf_POSTERIOR{i},0,1)/sqrt(M);
    pValue = zeros(1,X_num(i));
    for j = 1:X_num(i)
        [~,pValue(j)] = ttest(POSTERIOR{i}(:,j),shuf_POSTERIOR{i}(:,j));
    end
    hypothesis = fdrCorr(pValue,0.01);
    pos_sig = hypothesis & ave_POSTERIOR > ave_shuf_POSTERIOR;
    neg_sig = hypothesis & ave_POSTERIOR < ave_shuf_POSTERIOR;
    y = ones(1,X_num(i));
    scatter(X{i}(pos_sig),0.98*y(pos_sig),'s','filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    scatter(X{i}(neg_sig),0.02*y(neg_sig),'s','filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    patch([X{i} fliplr(X{i})],[ave_POSTERIOR+sem_POSTERIOR ...
        fliplr(ave_POSTERIOR-sem_POSTERIOR)],...
    'w','edgealpha',0.5); 
    patch([X{i} fliplr(X{i})],[ave_shuf_POSTERIOR+sem_shuf_POSTERIOR ...
        fliplr(ave_shuf_POSTERIOR-sem_shuf_POSTERIOR)],...
    'w','EdgeColor',[160 160 160]/255,'edgealpha',0.5);
    p1 = plot(X{i},ave_POSTERIOR,'Color','k','LineStyle','-','LineWidth',2);
    p2 = plot(X{i},ave_shuf_POSTERIOR,'Color','k','LineStyle','--','LineWidth',2);
    if i == 2
        legend([p1,p2],{'Attention cue loc','Control'},'Location','South');
        legend('boxoff');
    end
    xlim([tw{i}(1),tw{i}(end)]);
    ylim([0 1]);
    xlabel(timeLabel{i});
    ylabel('Posterior probability');
%     set(gca,'fontsize',14);
    hold on;
end
 
fig2 = figure('Position',[500 500 800 400]);
for i = 1:epoch_num    
    subplot(1,epoch_num,i);hold on;
    M= size(P,1);
    ave_POSTERIOR = mean(POSTERIOR_f{i},1);
    sem_POSTERIOR = std(POSTERIOR_f{i},0,1)/sqrt(M);
    ave_shuf_POSTERIOR = mean(shuf_POSTERIOR_f{i},1);
    sem_shuf_POSTERIOR = std(shuf_POSTERIOR_f{i},0,1)/sqrt(M);
    pValue = zeros(1,X_num(i));
    for j = 1:X_num(i)
        [~,pValue(j)] = ttest(POSTERIOR_f{i}(:,j),shuf_POSTERIOR_f{i}(:,j));
    end
    hypothesis = fdrCorr(pValue,0.01);
    pos_sig = hypothesis & ave_POSTERIOR > ave_shuf_POSTERIOR;
    neg_sig = hypothesis & ave_POSTERIOR < ave_shuf_POSTERIOR;
    y = ones(1,X_num(i));
    scatter(X{i}(pos_sig),0.98*y(pos_sig),'s','filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    scatter(X{i}(neg_sig),0.02*y(neg_sig),'s','filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    patch([X{i} fliplr(X{i})],[ave_POSTERIOR+sem_POSTERIOR ...
        fliplr(ave_POSTERIOR-sem_POSTERIOR)],...
    'w','edgealpha',0.5); 
    patch([X{i} fliplr(X{i})],[ave_shuf_POSTERIOR+sem_shuf_POSTERIOR ...
        fliplr(ave_shuf_POSTERIOR-sem_shuf_POSTERIOR)],...
    'w','EdgeColor',[160 160 160]/255,'edgealpha',0.5);
    p1 = plot(X{i},ave_POSTERIOR,'Color','k','LineStyle','-','LineWidth',2);
    p2 = plot(X{i},ave_shuf_POSTERIOR,'Color','k','LineStyle','--','LineWidth',2);
    if i == 2
        legend([p1,p2],{'Attention cue loc','Control'});
        legend('boxoff');
    end
    xlim([tw{i}(1),tw{i}(end)]);
    ylim([0 1]);
    xlabel(timeLabel{i});
    ylabel('Posterior probability');
%     set(gca,'fontsize',14);
    hold on;
end
end