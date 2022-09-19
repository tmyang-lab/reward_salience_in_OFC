function fig = plotEnsembleProbCrossTemp(ensembleNeuron)
%plotEnsembleProbCrossTemp Plot heatmap of posterior probability in 
% cross-temporal decoding
% Input: 
%    ensembleNeuron: is a structure, containing LDA results of neuronal
%    ensemble
% Output:
%    fig: heatmap

% Extract variables
P = ensembleNeuron.P_c;
shuf_P = ensembleNeuron.shuf_P_c;

% Plot figure
% Properties
% Bin width
bin_width = 25;
% Step 
bin_step = 10;
% Time epoch
timeLabel{1} = 'Time from stim on (ms)';
timeLabel{2} = 'Time from lum change (ms)'; 
tw{1} = [-400 600];
tw{2} = [-400 400];
% From bottom left to bottom right: train epoch 1 test epoch 1, train epoch
% 1 test epoch 2; from top left to top right: train epoch 1 test epoch 2,
% train epoch 2 test epoch 2
epoch_num = 2;
pos = [3 1;4 2];
% X is the time of each epoch align to imgage on and luminance change
% X_num is the length of each time epoch
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
% There are four matrix, each row is the training time, col is the testing
% time
POSTERIOR = cell(epoch_num,epoch_num);
shuf_POSTERIOR = cell(epoch_num,epoch_num);
for i = 1:epoch_num
    for j = 1:epoch_num
        POSTERIOR{i,j} = P(:,X_index{i},X_index{j});
        shuf_POSTERIOR{i,j} = shuf_P(:,X_index{i},X_index{j});
    end
end
fig = figure('Position',[100,100,900,900]);
for i = 1:epoch_num
    for j = 1:epoch_num
        subplot(epoch_num,epoch_num,pos(i,j)); hold on;                 
        ave_POSTERIOR = reshape(mean(POSTERIOR{i,j},1),[X_num(i),X_num(j)]);        
        image(X{i},X{j},ave_POSTERIOR','CDataMapping','scaled');
        caxis([0.2 0.8]);
        colorbar;
        colormap(viridis);
        xlim([tw{i}(1),tw{i}(end)]);
        ylim([tw{j}(1),tw{j}(end)]);
        xlabel(timeLabel{i});
        ylabel(timeLabel{j});
        set(gca,'FontSize',12);
        hold off;
    end
end

end