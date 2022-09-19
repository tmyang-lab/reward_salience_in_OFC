function fig = plotAccuracySubjD(x)
%plotAccuracySubjD   Plot accuracy in subject D.
% Input: x is a session number by 1 cell array, each cell is a trial
% number by variable number matrix. 
% Output: figure

m = length(x);
% Hit rate in valid-salient trial and valid-nonsalient trial
hr_v_s = zeros(m,1);
hr_v_ns = zeros(m,1);
% Hit rate in invalid-salient trial and invalid-nonsalient trial
hr_inv_s = zeros(m,1);
hr_inv_ns = zeros(m,1);
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
    % Miss trial
    miss_trial = x{i}(:,6);    
    % Hit rate in valid-salient trial and valid non-salient trial
    hr_v_s(i) = sum(dou_v_trial & s_trial & hit_trial)...
        /sum(dou_v_trial & s_trial & (hit_trial | miss_trial));
    hr_v_ns(i) = sum(dou_v_trial & ns_trial & hit_trial)...
        /sum(dou_v_trial & ns_trial & (hit_trial | miss_trial));
   % Hit rate in invalid-salient trial and invalid-nonsalient trial
    hr_inv_s(i) = sum(dou_inv_trial & s_trial & hit_trial)...
        /sum(dou_inv_trial & s_trial & (hit_trial | miss_trial));
    hr_inv_ns(i) = sum(dou_inv_trial & ns_trial & hit_trial)...
        /sum(dou_inv_trial & ns_trial & (hit_trial | miss_trial));
end

% Combine hit rate
accuracy_vInv_sNs(:,1,1) = hr_v_s;
accuracy_vInv_sNs(:,1,2) = hr_v_ns;
accuracy_vInv_sNs(:,2,1) = hr_inv_s;
accuracy_vInv_sNs(:,2,2) = hr_inv_ns;

% Plot figure
fig = figure; hold on;
X = 1:1:size(accuracy_vInv_sNs,2);
ave_accuracy_vInv_sNs = squeeze(mean(accuracy_vInv_sNs,1));
sem_accuracy_vInv_sNs = squeeze(std(accuracy_vInv_sNs,0,1))/sqrt(m);
% bar graph
b = bar(X,ave_accuracy_vInv_sNs,0.5,'FaceColor','w','EdgeColor','k','LineWidth',2);
drawnow;
% Error bar
x = zeros(2,2);
y = zeros(2,2);
for i = 1:2
    x(:,i) = bsxfun(@plus, b(1).XData, [b(i).XOffset]');
    y(:,i) = b(i).YData;
end
errorbar(x, y, sem_accuracy_vInv_sNs,'.','Color','k','LineWidth',2);
b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(1).CData = [135, 54, 0;135, 54, 0]/255;
b(2).CData = [237, 187, 153;237, 187, 153]/255;
legend('Salient','Non-salient');
legend('boxoff');
ylabel('Hit rate'); 
xticks(X);
xticklabels({'Valid','Invalid'});
ylim([0.5 1]); 
% set(gca,'FontSize',14);
hold off;

end