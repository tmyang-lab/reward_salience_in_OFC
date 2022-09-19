function fig = plotAccuracySubjG(x)
%plotAccuracySubjG   Plot accuracy in subject G.
% Input: x is a session number by 1 cell array, each cell is a trial
% number by variable number matrix. 
% Output: figure

m = length(x);
% Hit rate in target-salient trial and target-nonsalient trial
hr_t_s = zeros(m,1);
hr_t_ns = zeros(m,1);
% Correct rejection rate in distractor-only-salient trial and
% distractor-only-nonsalient trial or distractor+target-salient trial and
% distractor+target-nonsalient trial
cr_dodt_s = zeros(m,1);
cr_dodt_ns = zeros(m,1);
% Hit rate in distractor+target-salient trial and 
% distractor+target-nonsalient trial
hr_dt_s = zeros(m,1);
hr_dt_ns = zeros(m,1);
for i = 1:m
    % target-double-stimuli trial
    dou_t_trial = x{i}(:,1) == 2;
    % distractor-only-double-stimuli trial
    dou_do_trial = x{i}(:,1) == 3;
    % distractor+target-double-stimuli trial
    dou_dt_trial = x{i}(:,1) == 4;   
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
    s_trial = ((dou_t_trial | dou_dt_trial)  & (cue_value > uncue_value)) | (dou_do_trial & (cue_value < uncue_value));
    % trial with luminance change on non-salient stimulus
    ns_trial = ((dou_t_trial | dou_dt_trial) & (cue_value < uncue_value)) | (dou_do_trial & (cue_value > uncue_value));
    % hit trial
    hit_trial = x{i}(:,5);
    % miss trial
    miss_trial = x{i}(:,6); 
    % false alarm trial
    fa_trial = x{i}(:,7);
    % correct rejection trial
    cr_trial = x{i}(:,8);
    % Hit rate in target-salient trial and target-nonsalient trial
    hr_t_s(i) = sum(dou_t_trial & s_trial & hit_trial)...
        /sum(dou_t_trial & s_trial & (hit_trial | miss_trial));  
    hr_t_ns(i) = sum(dou_t_trial & ns_trial & hit_trial)...
        /sum(dou_t_trial & ns_trial & (hit_trial | miss_trial)); 
    % Correct rejection rate in distractor-only-salient trial and 
    % distractor-only-nonsalient trial or distractor+target-salient trial 
    % and distractor+target-nonsalient trial
    cr_dodt_s(i) = sum((dou_do_trial & s_trial & cr_trial) | (dou_dt_trial & ns_trial & cr_trial)) ...
        /sum((dou_do_trial & s_trial & (cr_trial | fa_trial)) | (dou_dt_trial & ns_trial & (cr_trial | fa_trial)));
    cr_dodt_ns(i) = sum((dou_do_trial & ns_trial & cr_trial) | (dou_dt_trial & s_trial & cr_trial)) ...
        /sum((dou_do_trial & ns_trial & (cr_trial | fa_trial)) | (dou_dt_trial & s_trial & (cr_trial | fa_trial)));
    % Hit rate in distractor+target-salient trial and
    % distractor+target-nonsalient trial
    hr_dt_s(i) = sum(dou_dt_trial & s_trial & hit_trial)...
        /sum(dou_dt_trial & s_trial & (hit_trial | miss_trial));   
    hr_dt_ns(i) = sum(dou_dt_trial & ns_trial & hit_trial)...
        /sum(dou_dt_trial & ns_trial & (hit_trial | miss_trial));
end

% Combine hit rate and correct rejection rate
accuracy_td_sNs(:,1,1) = hr_t_s;
accuracy_td_sNs(:,1,2) = hr_t_ns;
accuracy_td_sNs(:,2,1) = cr_dodt_s;
accuracy_td_sNs(:,2,2) = cr_dodt_ns;
accuracy_td_sNs(:,3,1) = hr_dt_s;
accuracy_td_sNs(:,3,2) = hr_dt_ns;

% Plot figure
fig = figure; hold on;
X = 1:1:size(accuracy_td_sNs,2);
ave_accuracy_td_sNs = squeeze(mean(accuracy_td_sNs,1));
sem_accuracy_td_sNs = squeeze(std(accuracy_td_sNs,0,1))/sqrt(m);
% bar graph
b = bar(X,ave_accuracy_td_sNs,0.5,'FaceColor','w','EdgeColor','k','LineWidth',2);
drawnow;
% Error bar
x = zeros(3,2);
y = zeros(3,2);
for i = 1:2
    x(:,i) = bsxfun(@plus, b(1).XData, [b(i).XOffset]');
    y(:,i) = b(i).YData;
end
errorbar(x, y, sem_accuracy_td_sNs,'.','Color','k','LineWidth',2);
b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(1).CData = [135, 54, 0;135, 54, 0;135, 54, 0]/255;
b(2).CData = [237, 187, 153;237, 187, 153;237, 187, 153]/255;
legend('Salient','Non-salient','location','North');
legend('boxoff');
ylabel('Hit/Corr rej rate'); 
xticks(X);
xticklabels({'Tgt(hit)','Dist/Dist+Tgt(CR)','Dist+Tgt(hit)'});
ylim([0.5 1]); 
% set(gca,'FontSize',14);
hold off;

end