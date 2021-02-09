% Script for a two part analysis. 

% 1) Analyse which reactions need alterations for the flux phenotype of EPI
% to take on a MES phenotype, representative of the changes in metabolism 
% that occur following EMT in breast.

% 2) Perform a flux enrichment analysis (using COBRA's FEA function) to 
% identify which reaction sets (within Recon 2.04) are significantly
% overrepresented in the list of reactions from step 1. 

% Sigurdur Karvelsson

%% Perform the EMT-linked reaction analysis in the proteomic GSMMs
initCobraToolbox;
load('sampled_proteomic_models.mat');
flux_E = FluxVectorize(modelENew,samplesE_32hours.points);
flux_M = FluxVectorize(modelMNew,samplesM_32hours.points);

% Check the altered reactions (MOMA-based test):
[relaxed_modelE,~,~,~,rxns_EMT] = relax_rxns(modelENew,[],[],[],flux_M,1,1);

% Check over-representation of reaction sets within the relaxed reactions:
id_rxns = find(ismember(modelENew.rxns,rxns_EMT(:,1)));
fea_res = FEA(modelENew,id_rxns,'subSystems');
T = cell2table(fea_res)
writetable(T,'EMT_reactions_enrichment_analysis.txt','Delimiter','\t')

