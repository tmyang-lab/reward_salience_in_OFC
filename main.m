% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Wenyi Zhang
% Date: 09/06/2022
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add to the path
mfPath = mfilename('fullpath');
addpath(genpath(fileparts(mfPath)));

%% Load data
% Behavioral data of monkey D
bhvD = load('data/bhvD.mat');
% Behavioral data of monkey G
bhvG = load('data/bhvG.mat');
% Example OFC neuron
exampleOFC = load('data/exampleOFC.mat');
% Exampe DLPFC neuron
exampleDLPFC = load('data/exampleDLPFC.mat');
% Population OFC neurons
populationOFC = load('data/populationOFC.mat');
% Population DLPFC neurons
populationDLPFC = load('data/populationDLPFC.mat');
% DLPFC neuronal ensemble 
ensembleDLPFC = load('data/ensembleDLPFC.mat');
% OFC neuronal ensemble 
ensembleOFC = load('data/ensembleOFC.mat');

%% Figure 1 b, c, e, f
fig1b = plotAccuracySubjD(bhvD.task_info);
fig1c = plotRTSubjD(bhvD.task_info);
fig1e = plotAccuracySubjG(bhvG.task_info);
fig1f = plotRTSubjG(bhvG.task_info);

%% Figure 2 a£¬b, c, d, e
fig2a = plotExampleDLPFCResp(exampleDLPFC);
% Whether recalculate the LDA results
isRecalLDA = false;
if isRecalLDA
    % Posterior of attention cue location in cross-temporal decoding
    [ensembleDLPFC.P_c,ensembleDLPFC.shuf_P_c] ...
        = LDARes(ensembleDLPFC.norm_fr,...
        ensembleDLPFC.att_cue_loc,ensembleDLPFC.CV,ensembleDLPFC.UCV,1);
    % Posterior of attention cue location trained and tested with same time bin
    [ensembleDLPFC.P,ensembleDLPFC.shuf_P] = LDARes(ensembleDLPFC.norm_fr,...
        ensembleDLPFC.att_cue_loc,ensembleDLPFC.CV,ensembleDLPFC.UCV,2);
    % Posterior of attention cue location trained before stimulus onset
    [ensembleDLPFC.P_f,ensembleDLPFC.shuf_P_f] = LDARes(ensembleDLPFC.norm_fr,...
        ensembleDLPFC.att_cue_loc,ensembleDLPFC.CV,ensembleDLPFC.UCV,3);
    % Posterior of attention cue location in quick and slow trials
    [ensembleDLPFC.P_DLPFC_q,~] = LDARes(ensembleDLPFC.norm_fr_q,...
        ensembleDLPFC.att_cue_loc_qs,[],[],2,false);
    [ensembleDLPFC.P_DLPFC_s,~] = LDARes(ensembleDLPFC.norm_fr_s,...
        ensembleDLPFC.att_cue_loc_qs,[],[],2,false);
end
fig2b = plotEnsembleProbCrossTemp(ensembleDLPFC);
[fig2c,fig2d] = plotEnsembleProb(ensembleDLPFC);
fig2e = plotEnsembleProbQuickSlow(ensembleDLPFC);

%% Figure 3
[fig3a,fig3d,fig3g,fig3j] = plotExampleOFCResp(exampleOFC);
[fig3b,fig3e,fig3h,fig3k] = plotValueTunedResp(populationOFC.singleNeuronFrVar(populationOFC.posIdx));
[fig3c,fig3f,fig3i,fig3l] = plotValueTunedResp(populationOFC.singleNeuronFrVar(populationOFC.negIdx));

%% Figure 4
% Whether recalculate the linear regression results
isRecalLinear = false;
if isRecalLinear
    populationOFC.singleNeuronRegResults = regMultiVar(populationOFC.singleNeuronFrVar);  
end
[fig4a,fig4d] = plotCPDCrossTime(populationOFC.singleNeuronRegResults);
[fig4b,fig4e] = plotCPDScatter(populationOFC.singleNeuronRegResults);
[fig4c,fig4f] = plotCPDHistogrm(populationOFC.singleNeuronRegResults);

%% Figure 5
fig5a = plotValueTunedPaired8Resp(populationOFC.singleNeuronFrVar(populationOFC.posIdx));
fig5b = plotValueTunedPaired8Resp(populationOFC.singleNeuronFrVar(populationOFC.negIdx));

%% Figure 6
% isRecalLinear = false;
isRecalLinear = true;
if isRecalLinear    
    populationDLPFC.singleNeuronRegResults = regMultiVar(populationDLPFC.singleNeuronFrVar);
end
[fig6a,fig6d] = plotCPDCrossTime(populationDLPFC.singleNeuronRegResults);
[fig6b,fig6e] = plotCPDScatter(populationDLPFC.singleNeuronRegResults);
[fig6c,fig6f] = plotCPDHistogrm(populationDLPFC.singleNeuronRegResults);

%% Figure 7
regAttCueSV_OFC = populationOFC.singleNeuronRegResults.aveCueOff_attCue_SV_NSV;
regAttCueSV_DLPFC = populationDLPFC.singleNeuronRegResults.aveCueOff_attCue_SV_NSV;
[fig7,fig7_insert] = plotCPD2Areas(regAttCueSV_OFC,regAttCueSV_DLPFC);
