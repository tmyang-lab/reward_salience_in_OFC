function fig = plotRTSubjD(x)
%plotRTSubjD   Plot reaction time in subject D.
% Input: x is a session number by 1 cell array, each cell is a trial
% number by variable number matrix. 
% Output: figure

m = length(x);
% Reaction time in valid-salient trial and valid-nonsalient trial
rt_v_s = zeros(m,1);
rt_v_ns = zeros(m,1);
% Reaction time in invalid-salient trial and invalid-nonsalient trial
rt_inv_s = zeros(m,1);
rt_inv_ns = zeros(m,1);
for i = 1:m
    % valid-double-stimuli trial
    dou_v_trial = x{i}(:,1) == 2;
    % invalid-double-stimuli trial
    dou_inv_trial = x{i}(:,1) == 3;  
    % attention cue location
    att_cue_loc = x{i}(:,2);
    % left value
    l_value = x{i}(:,3);
    % right value
    r_value = x{i}(:,4);
    % cued value
    cue_value = (att_cue_loc == -1).*r_value + (att_cue_loc == 1).*l_value;
    % uncued value
    uncue_value = (att_cue_loc == -1).*l_value + (att_cue_loc == 1).*r_value;
    % trial with luminance change on salient stimulus
    s_trial = (dou_v_trial & (cue_value > uncue_value)) | (dou_inv_trial & (cue_value < uncue_value));
    % trial with luminance change on non-salient stimulus
    ns_trial = (dou_v_trial & (cue_value < uncue_value)) | (dou_inv_trial & (cue_value > uncue_value));
    % Hit trial
    hit_trial = x{i}(:,5);   
    % Reaction time
    RT = x{i}(:,9);
    % Reaction time in valid-salient trial and valid non-salient trial
    rt_v_s(i) = mean(RT(dou_v_trial & s_trial & hit_trial));
    rt_v_ns(i) = mean(RT(dou_v_trial & ns_trial & hit_trial));
    % Reaction time in invalid-salient trial and invalid-nonsalient trial
    rt_inv_s(i) = mean(RT(dou_inv_trial & s_trial & hit_trial));
    rt_inv_ns(i) = mean(RT(dou_inv_trial & ns_trial & hit_trial));
end

% Combine reaction time in valid-salient trial, valid-nonsalient, 
% invalid-salient trial, and invalid-nonsalient trial
rt_vInv_sNs(:,1,1) = rt_v_s;
rt_vInv_sNs(:,1,2) = rt_v_ns;
rt_vInv_sNs(:,2,1) = rt_inv_s;
rt_vInv_sNs(:,2,2) = rt_inv_ns;

%% Plot figure
fig = figure; hold on;
X = 1:1:size(rt_vInv_sNs,2);
ave_rt_vInv_sNs = squeeze(mean(rt_vInv_sNs,1));
sem_rt_vInv_sNs = squeeze(std(rt_vInv_sNs,0,1))/sqrt(m);
% bar graph
b = bar(X,ave_rt_vInv_sNs,'FaceColor','w','EdgeColor','k','LineWidth',2);
drawnow;
% Error bar
x = zeros(2,2);
y = zeros(2,2);
for i = 1:2
    x(:,i) = bsxfun(@plus, b(1).XData, [b(i).XOffset]');
    y(:,i) = b(i).YData;
end
errorbar(x, y, sem_rt_vInv_sNs,'.','Color','k','LineWidth',2);
b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(1).CData = [135, 54, 0;135, 54, 0]/255;
b(2).CData = [237, 187, 153;237, 187, 153]/255;
legend('Salient','Non-salient','location','NorthWest');
legend('boxoff');
ylabel('Reaction time (ms)'); 
xticks(X);
xticklabels({'Valid','Invalid'});
ylim([150 250]); 
% set(gca,'FontSize',14);
hold off;

end