function [cP,shuf_cP] = LDARes(y, x, x1, x2, cross_temp, fit_via_res, train_period, varargin)
%LDARes Perform LDA on the residuals
% Perform LDA on y to predict x after excluding the influence of x1 and x2
% Inputs:
%   - (1) y: y is an M by T by N matrix, M is trial number, T is time
%            number, and N is neuron number.
%   - (2) x: x is an M by 1 vector. x is the variable predicted by the
%            model.
%   - (3) x1: x1 is an M by 1 vector. x1 is the variable that might
%         correlate with x.
%   - (4) x2: x2 is an M by 1 vector. x2 might also correlate with x.
%   - (5) cross_temp: cross_temp could be 1, 2, or 3. 1 means 
%         cross-temporal decoding. 2 means training and testing in same
%         time bin. 3 means training with train_period responses and
%         testing with all period responses
%   - (6) fit_via_res (default: true): if true, LDA performed on the
%         residuals after ANOVA with x1 and x2
%   - (7) train_period (optional): if cross_temp is 3, trained with
%         responses during train_period
% Outputs:
%   - (1) cP: cP is an M by T matrix. Posterior probability in true data
%   - (2) shuf_cP: shuf_cP is an M by T matrix. Posterior probability in 
%         shuffled data

if nargin < 5 || isempty(intersect(cross_temp,[1,2,3])), cP = []; shuf_cP = []; end
if nargin < 6 || isempty(fit_via_res), fit_via_res = true; end
if nargin < 7 || isempty(train_period), train_period = [-200,-50]; end

% Time of each data points
bin_width = 25;
t_step = 10;
tw = [-400 600];
X_ini = 1+tw(1):t_step:tw(end)-bin_width+1;
X_end = tw(1)+bin_width:t_step:tw(end);
X_idx = (X_ini + X_end)/2;

% Fit model vis residuals
M = size(y,1);
T = size(y,2);
N = size(y,3);
if fit_via_res
    y_fit = zeros(M,T,N);
    for i = 1:T
        fprintf(['>>>> Calculating residuals of time bin %d (',num2str(T),')... \n'],i);
        for j = 1:N
            [~,~,stats] = anovan(y(:,i,j),{x1(:,j),x2(:,j)},'display','off');
            y_fit(:,i,j) = stats.resid;
        end
    end
else
    y_fit = y;
end

% LOOV
train_idx = zeros(M,M-1);
test_idx = zeros(M,1);
for k = 1:M
    test_idx(k,:) = k;
    train_idx(k,:) = setdiff(1:1:M,test_idx(k));
end

% LDA
if cross_temp == 1  % Cross-temporal decoding 
    fprintf('>>>> Performing cross-temporal decoding \n');
    cP = zeros(M,T,T);
    shuf_cP = zeros(M,T,T);
    % perform pca
    yfit_reshape = reshape(y_fit,[M*T,N]);
    [~,score,~,~,explained,~] = pca(yfit_reshape);
    cum_explained = cumsum(explained/sum(explained));
    pc_num = sum(cum_explained < 0.7) + 1;
    yfit_pca = reshape(score(:,1:pc_num),[M,T,pc_num]);
    % perform lda    
    for i = 1:M
        fprintf(['>>>> Calculating posterior of trial %d (',num2str(M),')... \n'],i);
        training_label = x(train_idx(i,:));
        testing_label = x(test_idx(i,:));
        shuf_x = x(randperm(M));
        shuf_testing_label = shuf_x(test_idx(i,:));
        for j = 1:T
            for k = 1:T
                training_data = reshape(yfit_pca(train_idx(i,:),j,:),[M-1,pc_num]);
                testing_data = reshape(yfit_pca(test_idx(i,:),k,:),[1,pc_num]);
                [~,~,POSTERIOR] = classify(testing_data,training_data,training_label);
                cP(i,j,k) = POSTERIOR([-1 1] == testing_label);
                shuf_cP(i,j,k) = POSTERIOR([-1 1]  == shuf_testing_label);
            end
        end
    end
elseif cross_temp == 2 % Trained and tested with same time bin  
    fprintf('>>>> Training and testing the decoder with reponses in the same time bin \n');
    cP = zeros(M,T);
    shuf_cP = zeros(M,T);
    for i = 1:M 
        fprintf(['>>>> Calculating posterior of trial %d (',num2str(M),')... \n'],i);
        training_label = x(train_idx(i,:));
        testing_label = x(test_idx(i,:));
        shuf_x = x(randperm(M));
        shuf_testing_label = shuf_x(test_idx(i,:));
        for j = 1:T
            % perform pca
            yfit_reshape = reshape(y_fit(:,j,:),[M,N]);
            [~,score,~,~,explained,~] = pca(yfit_reshape);
            cum_explained = cumsum(explained/sum(explained));
            pc_num = sum(cum_explained < 0.7) + 1;
            yfit_pca = score(:,1:pc_num);
            training_data = yfit_pca(train_idx(i,:),:);  
            testing_data = yfit_pca(test_idx(i,:),:);
            [~,~,POSTERIOR] = classify(testing_data,training_data,training_label);
            cP(i,j) = POSTERIOR([-1 1] == testing_label);
            shuf_cP(i,j) = POSTERIOR([-1 1] == shuf_testing_label);
        end
    end   
elseif cross_temp == 3 % Trained before stimulus onset
    fprintf('>>>> Training the decoder with reponses before stimulus onset \n');
    cP = zeros(M,T);
    shuf_cP = zeros(M,T); 
    % Perform pca
    yfit_reshape = reshape(y_fit,[M*T,N]);
    [~,score,~,~,explained,~] = pca(yfit_reshape);
    cum_explained = cumsum(explained/sum(explained));
    pc_num = sum(cum_explained < 0.7) + 1;
    yfit_pca = reshape(score(:,1:pc_num),[M,T,pc_num]);
    for i = 1:M 
        fprintf(['>>>> Calculating posterior of trial %d (',num2str(M),')... \n'],i);
        training_label = x(train_idx(i,:));
        testing_label = x(test_idx(i,:));
        shuf_x = x(randperm(M));
        shuf_testing_label = shuf_x(test_idx(i,:));        
        training_data = reshape(mean(yfit_pca(train_idx(i,:),...
            X_idx>=train_period(1) & X_idx<train_period(end),:),2),[M-1,pc_num]);        
        for j = 1:T
            % Perform lda
            testing_data = reshape(yfit_pca(test_idx(i,:),j,:),[1,pc_num]);
            [~,~,POSTERIOR] = classify(testing_data,training_data,training_label);
            cP(i,j) = POSTERIOR([-1 1] == testing_label);
            shuf_cP(i,j) = POSTERIOR([-1 1] == shuf_testing_label);
        end    
    end  
end
end