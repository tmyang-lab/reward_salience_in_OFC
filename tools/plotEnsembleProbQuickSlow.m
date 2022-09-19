function fig = plotEnsembleProbQuickSlow(ensembleNeuron)
%plotEnsembleProbQuickSlow Plot shaded area of posterior probability in
% quickly and slowly reacted trials

% Input: 
%    ensembleNeuron: is a structure, containing LDA results of neuronal
%    ensemble
% Output:
%    fig

% Extract variables
P_q = ensembleNeuron.P_q;
P_s = ensembleNeuron.P_s;
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
POSTERIOR_q = cell(1,epoch_num);
POSTERIOR_s = cell(1,epoch_num);
for i = 1:epoch_num
        POSTERIOR_q{i} = P_q(:,X_index{i});
        POSTERIOR_s{i} = P_s(:,X_index{i});
end

% Plot figure
fig = figure('Position',[500 500 800 400]);
for i = 1:epoch_num    
    subplot(1,epoch_num,i);hold on;
    M= size(P_q,1);
    ave_POSTERIOR_q = mean(POSTERIOR_q{i},1);
    sem_POSTERIOR_q = std(POSTERIOR_q{i},0,1)/sqrt(M);
    ave_POSTERIOR_s = mean(POSTERIOR_s{i},1);
    sem_POSTERIOR_s = std(POSTERIOR_s{i},0,1)/sqrt(M);
    pValue = zeros(1,X_num(i));
    for j = 1:X_num(i)
        [~,pValue(j)] = ttest(POSTERIOR_q{i}(:,j),POSTERIOR_s{i}(:,j));
    end
    hypothesis = fdrCorr(pValue,0.01);
    pos_sig = hypothesis & ave_POSTERIOR_q > ave_POSTERIOR_s;
    neg_sig = hypothesis & ave_POSTERIOR_q < ave_POSTERIOR_s;
    y = ones(1,X_num(i));
    scatter(X{i}(pos_sig),0.98*y(pos_sig),'s','filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    scatter(X{i}(neg_sig),0.02*y(neg_sig),'s','filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    patch([X{i} fliplr(X{i})],[ave_POSTERIOR_q+sem_POSTERIOR_q ...
        fliplr(ave_POSTERIOR_q-sem_POSTERIOR_q)],...
    'w','edgealpha',0.5); 
    patch([X{i} fliplr(X{i})],[ave_POSTERIOR_s+sem_POSTERIOR_s ...
        fliplr(ave_POSTERIOR_s-sem_POSTERIOR_s)],...
    'w','EdgeColor',[160 160 160]/255,'edgealpha',0.5);
    p1 = plot(X{i},ave_POSTERIOR_q,'Color','k','LineStyle','-','LineWidth',2);
    p2 = plot(X{i},ave_POSTERIOR_s,'Color',[128 128 128]/255,'LineStyle','-','LineWidth',2);
    if i == 2
        legend([p1,p2],{'Quick','Slow'},'Location','South');
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