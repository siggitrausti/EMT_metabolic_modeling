% Script to create the different data-type constrained GSMMs from
% iBreast2886. This is a modified version of the code published in
% Halldorsson et al., 2017 (doi: 10.1016/j.canlet.2017.03.019). 

% Sigurdur Karvelsson

%% Load the iBreast2886 model (from Halldorsson et al.)
initCobraToolbox;
load('EMTmodel.mat');
steps = 5800 % steps for each point. 
time = 32*60*60; % 32 hours in seconds
model.csense(1:size(model.S,1),1) = 'E'; % Lacking the constraint sense.
% Use equality, which is the default
[sampleMetaOut, MF] = gpSampler(model,steps,[],time);
save('sampled_iBreast2886.mat','sampleMetaOut');

%% Load necessary measurements
% From Halldorsson et al.:
load('EMTmodel.mat');
load('MetabolomicsFluxEpi0to48.mat');
load('MetabolomicsFluxMes0to48.mat');
load('Met2EXrxns.mat');

% Other:
load('ID_conversion_gene_names_to_ENTREZ.mat');
load('sampled_iBreast2886.mat');

%% 1) Creation of the only media-constrained iBreast2886 GSMMs
% The amount of flux downregulation (1 corresponds to knockout)
tol = 1e-5; 
downRegFrac=0.7;
downRegFrac_abs=0.7;

nrxns=length(model.rxns);
newlb=zeros(nrxns,1);
newub=zeros(nrxns,1);
lower=zeros(nrxns,1);
upper=zeros(nrxns,1);

for i=1:length(model.rxns)
    % Default values correspond to the original upper/lower bounds
    newlb(i)=model.lb(i);
    newub(i)=model.ub(i);
    if newlb(i) == newub(i)
        % Fixed reaction, possibly disabled
        fprintf('%d\t%s has fixed upper/lower bounds (%1.4f)\n', i, model.rxns{i}, model.lb(i))
        continue
    end
    flux=sampleMetaOut.points(i,:);
    lower(i)=min(flux);
    upper(i)=max(flux);
    newlb(i)=min(0,prctile(flux(flux<tol), downRegFrac*100));
    newub(i)=max(0,prctile(flux(flux>tol), (1 - downRegFrac)*100));
end

% Create the EPI and MES models
modelNew = model
modelNew.lb=lower;
modelNew.ub=upper;

% Set all bounds to MS data
modelE=modelNew;
modelM=modelNew;

forcedExchangesE = [];
forcedExchangesM = [];
percentError = 0.1;
Ex = [];

for i = 2:size(fluxdataEpi0to48,1)
    Ex = MetNames2ExRxns((find(ismember(MetNames2ExRxns(:,1),fluxdataEpi0to48(i,1)))),2);
    if ~isempty(find(ismember(modelE.rxns,Ex)))
        
        forcedExchangesE = cat(1,forcedExchangesE,[Ex,cell2mat(fluxdataEpi0to48(i,2))-(abs(cell2mat(fluxdataEpi0to48(i,2))*percentError)),cell2mat(fluxdataEpi0to48(i,2))+(abs(cell2mat(fluxdataEpi0to48(i,2))*percentError))]);
        modelE=changeRxnBounds(modelE, Ex, cell2mat(fluxdataEpi0to48(i,2))-abs((cell2mat(fluxdataEpi0to48(i,2))*percentError)), 'l');
        modelE=changeRxnBounds(modelE, Ex, cell2mat(fluxdataEpi0to48(i,2))+abs((cell2mat(fluxdataEpi0to48(i,2))*percentError)), 'u');
       
        forcedExchangesM = cat(1,forcedExchangesM,[Ex,cell2mat(fluxdataMes0to48(i,2))-abs((cell2mat(fluxdataMes0to48(i,2))*percentError)),cell2mat(fluxdataMes0to48(i,2))+abs((cell2mat(fluxdataMes0to48(i,2))*percentError))]);
        modelM=changeRxnBounds(modelM, Ex, cell2mat(fluxdataMes0to48(i,2))-abs((cell2mat(fluxdataMes0to48(i,2))*percentError)), 'l');
        modelM=changeRxnBounds(modelM, Ex, cell2mat(fluxdataMes0to48(i,2))+abs((cell2mat(fluxdataMes0to48(i,2))*percentError)), 'u');

    else
    end
end

% To check for the metabolites which have huge errors in MS data
% Glutamine measured in a different experiment for accuracy
modelE.lb(find(ismember(modelE.rxns,'EX_pyr(e)'))) = -10;
modelE.ub(find(ismember(modelE.rxns,'EX_pyr(e)'))) = 0;
modelE.lb(find(ismember(modelE.rxns,'EX_lac_D(e)'))) = 0;
modelE.ub(find(ismember(modelE.rxns,'EX_lac_D(e)'))) = 100;

modelM.lb(find(ismember(modelM.rxns,'EX_pyr(e)'))) = -10;
modelM.ub(find(ismember(modelM.rxns,'EX_pyr(e)'))) = 0;
modelM.lb(find(ismember(modelM.rxns,'EX_lac_D(e)'))) = 0;
modelM.ub(find(ismember(modelM.rxns,'EX_lac_D(e)'))) = 100;
idx_relax =find(ismember(model.rxns,model.rxns));

% Relax models to achieve feasibility:
FBAE = optimizeCbModel(modelE);
FBAE.f
FBAM = optimizeCbModel(modelM);
FBAM.f
if FBAE.f == 0 || FBAM.f == 0
    try
    [modelENew,r,d,f,rxns_relaxed] = relax_rxns(modelE,idx_relax,0.01,[],[],0);
    [modelMNew,r,d,f,rxns_relaxed] = relax_rxns(modelM,idx_relax,0.01,[],[],0);
    catch
    end
else
    modelENew = modelE
    modelMNew = modelM
end

% ----------------- FBA -------------------------
% Check for feasibility of relaxed models:
FBAE = optimizeCbModel(modelENew);
if FBAE.f < tol
    fprintf('*************************************************\n')
    fprintf('Infeasible EPI model!\n')
else
    fprintf('*************************************************\n')
    fprintf('Feasible EPI model achieved!\n')
end

FBAM = optimizeCbModel(modelMNew);
if FBAM.f < tol
    fprintf('Infeasible MES model!\n')
else
    fprintf('Feasible MES model achieved!\n')
end
% ----------------- GEA -------------------------
[grRatioE,~,~,~,hasEffectE] = singleGeneDeletion(modelENew,'FBA',modelENew.genes,'true');
lethalE = find(isnan(grRatioE));
[grRatioM,~,~,~,hasEffectM] = singleGeneDeletion(modelMNew,'FBA',modelMNew.genes,'true');
lethalM = find(isnan(grRatioM));
lethal_onlyM = model.genes(setdiff(lethalM,lethalE));
lethal_onlyE = model.genes(setdiff(lethalE,lethalM));

save('media_EMT_GSMMs.mat')
clear
%% 2) Creation of the microarray-constrained iBreast2886 GSMMs
% Load necessary measurements
% From Halldorsson et al.:
load('EMTmodel.mat');
load('MetabolomicsFluxEpi0to48.mat');
load('MetabolomicsFluxMes0to48.mat');
load('Met2EXrxns.mat');

% Other:
load('ID_conversion_gene_names_to_ENTREZ.mat');
load('sampled_iBreast2886.mat');

% Microarray (directly from Halldorsson et al.):
% Converted ILMN IDs to Entrez gene IDs using online converter DAVID

[NUMcnvData,TXTcnvData,RAWcnvData] = xlsread('converted_EntrezID_DAVID.xls', 1, 'A2:B14768');
[NUMdata,TXTdata,RAWdata] = xlsread('E_vs_M.xls',1,'A2:F13106');
Cut_off = 5;
[MicorArray_Recon2] = ExpData_to_RxnNames_STK(RAWdata,model,RAWcnvData,'recon2');
down_regulated = MicorArray_Recon2(find(cell2mat(MicorArray_Recon2(:,3)) < -(Cut_off)),:);
fprintf('Number of Recon2 Genes with decreased expression level: \t %d\n', length(unique(down_regulated(:,7))))
up_regulated = MicorArray_Recon2(find(cell2mat(MicorArray_Recon2(:,3)) > Cut_off),:);
fprintf('Number of Recon2 Genes with increased expression level: \t %d\n', length(unique(up_regulated(:,7))))

% Create the models:
tol = 1e-5;

% The amount of flux downregulation (1 corresponds to knockout)
downRegFrac=0.7;
downRegFrac_abs=0.7;

nrxns=length(model.rxns);
newlb=zeros(nrxns,1);
newub=zeros(nrxns,1);
lower=zeros(nrxns,1);
upper=zeros(nrxns,1);

for i=1:length(model.rxns)
    % Default values correspond to the original upper/lower bounds
    newlb(i)=model.lb(i);
    newub(i)=model.ub(i);
    if newlb(i) == newub(i)
        % Fixed reaction, possibly disabled
        fprintf('%d\t%s has fixed upper/lower bounds (%1.4f)\n', i, model.rxns{i}, model.lb(i))
        continue
    end
    flux=sampleMetaOut.points(i,:);
    lower(i)=min(flux);
    upper(i)=max(flux);
    newlb(i)=min(0,prctile(flux(flux<tol), downRegFrac*100));
    newub(i)=max(0,prctile(flux(flux>tol), (1 - downRegFrac)*100));
end

% Create the EPI and MES models
modelNew = model;
modelNew.lb=lower;
modelNew.ub=upper;

% Constrain using microarray data:
idx_M=affectedRxns(model,down_regulated(:,7));
idx_E=affectedRxns(model,up_regulated(:,7));

fprintf('EPI model: Number of affected reactions: %d\n', length(idx_E))
modelE=modelNew;
for i=1:length(idx_E)
    rxnE=idx_E(i);
    modelE=changeRxnBounds(modelE, modelE.rxns{rxnE}, newlb(rxnE),'l');
    modelE=changeRxnBounds(modelE, modelE.rxns{rxnE}, newub(rxnE),'u');
    fprintf('EPI\t%s\t%1.4f\t%1.4f\n', modelE.rxns{rxnE}, modelE.lb(rxnE), modelE.ub(rxnE))
end

for i=1:length(modelE.rxns)
    % Default values correspond to the original upper/lower bounds
    flux_E=sampleMetaOut.points(i,:);
    newlb_E(i)=min(0,prctile(flux_E(flux_E<tol), downRegFrac*100));
    newub_E(i)=max(0,prctile(flux_E(flux_E>tol), (1 - downRegFrac)*100));
end


fprintf('*************************************************\n')
fprintf('\nMES model: Number of affected reactions: %d\n', length(idx_M))
modelM=modelNew;
for i=1:length(idx_M)
    rxnM=idx_M(i);
    modelM=changeRxnBounds(modelM, modelM.rxns{rxnM}, newlb(rxnM),'l');
    modelM=changeRxnBounds(modelM, modelM.rxns{rxnM}, newub(rxnM),'u');
    fprintf('MES\t%s\t%1.4f\t%1.4f\n', modelM.rxns{rxnM}, modelM.lb(rxnM), modelM.ub(rxnM))
end

for i=1:length(modelM.rxns)
    % Default values correspond to the original upper/lower bounds
    flux_M = sampleMetaOut.points(i,:);
    newlb_M(i)=min(0,prctile(flux_M(flux_M<tol), downRegFrac*100));
    newub_M(i)=max(0,prctile(flux_M(flux_M>tol), (1 - downRegFrac)*100));
end

% Set all bounds to MS data
forcedExchangesE = [];
forcedExchangesM = [];
percentError = 0.1;
Ex = [];

for i = 2:size(fluxdataEpi0to48,1)
    Ex = MetNames2ExRxns((find(ismember(MetNames2ExRxns(:,1),fluxdataEpi0to48(i,1)))),2);
    if ~isempty(find(ismember(modelE.rxns,Ex)))
        
        forcedExchangesE = cat(1,forcedExchangesE,[Ex,cell2mat(fluxdataEpi0to48(i,2))-(abs(cell2mat(fluxdataEpi0to48(i,2))*percentError)),cell2mat(fluxdataEpi0to48(i,2))+(abs(cell2mat(fluxdataEpi0to48(i,2))*percentError))]);
        modelE=changeRxnBounds(modelE, Ex, cell2mat(fluxdataEpi0to48(i,2))-abs((cell2mat(fluxdataEpi0to48(i,2))*percentError)), 'l');
        modelE=changeRxnBounds(modelE, Ex, cell2mat(fluxdataEpi0to48(i,2))+abs((cell2mat(fluxdataEpi0to48(i,2))*percentError)), 'u');
       
        forcedExchangesM = cat(1,forcedExchangesM,[Ex,cell2mat(fluxdataMes0to48(i,2))-abs((cell2mat(fluxdataMes0to48(i,2))*percentError)),cell2mat(fluxdataMes0to48(i,2))+abs((cell2mat(fluxdataMes0to48(i,2))*percentError))]);
        modelM=changeRxnBounds(modelM, Ex, cell2mat(fluxdataMes0to48(i,2))-abs((cell2mat(fluxdataMes0to48(i,2))*percentError)), 'l');
        modelM=changeRxnBounds(modelM, Ex, cell2mat(fluxdataMes0to48(i,2))+abs((cell2mat(fluxdataMes0to48(i,2))*percentError)), 'u');

    else
    end
end

% To check for the metabolites which have huge errors in MS data
% Glutamine measured in a different experiment for accuracy
modelE.lb(find(ismember(modelE.rxns,'EX_pyr(e)'))) = -10;
modelE.ub(find(ismember(modelE.rxns,'EX_pyr(e)'))) = 0;
modelE.lb(find(ismember(modelE.rxns,'EX_lac_D(e)'))) = 0;
modelE.ub(find(ismember(modelE.rxns,'EX_lac_D(e)'))) = 100;

modelM.lb(find(ismember(modelM.rxns,'EX_pyr(e)'))) = -10;
modelM.ub(find(ismember(modelM.rxns,'EX_pyr(e)'))) = 0;
modelM.lb(find(ismember(modelM.rxns,'EX_lac_D(e)'))) = 0;
modelM.ub(find(ismember(modelM.rxns,'EX_lac_D(e)'))) = 100;
idx_relax =find(ismember(model.rxns,model.rxns));

% Relax models to achieve feasibility:
FBAE = optimizeCbModel(modelE);
FBAE.f
FBAM = optimizeCbModel(modelM);
FBAM.f
if FBAE.f == 0 || FBAM.f == 0
    try
    [modelENew,r,d,f,rxns_relaxed] = relax_rxns(modelE,idx_relax,0.01,[],[],0);
    [modelMNew,r,d,f,rxns_relaxed] = relax_rxns(modelM,idx_relax,0.01,[],[],0);
    catch
    end
else
    modelENew = modelE
    modelMNew = modelM
end

% ----------------- FBA -------------------------
% Check for feasibility of relaxed models:
FBAE = optimizeCbModel(modelENew);
if FBAE.f < tol
    fprintf('*************************************************\n')
    fprintf('Infeasible EPI model!\n')
else
    fprintf('*************************************************\n')
    fprintf('Feasible EPI model achieved!\n')
end

FBAM = optimizeCbModel(modelMNew);
if FBAM.f < tol
    fprintf('Infeasible MES model!\n')
else
    fprintf('Feasible MES model achieved!\n')
end

% ----------------- GEA -------------------------
[grRatioE,~,~,~,hasEffectE] = singleGeneDeletion(modelENew,'FBA',modelENew.genes,'true');
lethalE = find(isnan(grRatioE));
[grRatioM,~,~,~,hasEffectM] = singleGeneDeletion(modelMNew,'FBA',modelMNew.genes,'true');
lethalM = find(isnan(grRatioM));
lethal_onlyM = model.genes(setdiff(lethalM,lethalE));
lethal_onlyE = model.genes(setdiff(lethalE,lethalM));

save('microarray_EMT_GSMMs.mat')
clear
%% 3) Creation of the proteomic-constrained iBreast2886 GSMMs
% Load necessary measurements
% From Halldorsson et al.:
load('EMTmodel.mat');
load('MetabolomicsFluxEpi0to48.mat');
load('MetabolomicsFluxMes0to48.mat');
load('Met2EXrxns.mat');

% Other:
load('ID_conversion_gene_names_to_ENTREZ.mat');
load('sampled_iBreast2886.mat');

[~,~,raw] = xlsread('Analysis Results Qiong.xlsx') % proteomic data
[~,~,rawID] = xlsread('IDconversion_proteomic.xlsx') % id conversion file
raw2 = raw(:,[1,5,13,14,10,11,16]) 
raw3 = raw2;
% Remove any NaN values:
raw2(any(cellfun(@(x) any(isnan(x)),raw2),2),:) = [];
% Change zeros to Nan:
for i=1:(size(raw2,1)*size(raw2,2))
    if cell2mat(raw2(i)) == 0;
        raw2(i) = {NaN};
    else
        raw2(i) = raw2(i);
    end
end
% Log transform data (to get normailzation of data)
log_data = log(cell2mat(raw2(2:end,2:7)))
data_transformed = [raw2(2:end,1) num2cell(log_data)]

% Perform a student's t-test:
pvalues = mattest(cell2mat(data_transformed(:,2:4)), cell2mat(data_transformed(:,5:7)));

% Correct for multiple comparisons:
[FDR, Q] = mafdr(pvalues);
FDR_test = linspace(0.01,0.1,10)'

% This is from another script for sensitivity analysis using different
% significance threshold values. Don't get confused by this!
prot_sens = num2cell(zeros(length(FDR_test),3))
for i=1:10
    prot_sens{i,2} = data_transformed(find(FDR < FDR_test(i)),1);
end
prot_sens(:,1) = num2cell(FDR_test(:))

% Now, to map IDs to entrez:
[Expression_Model] = ExpData_to_RxnNames_STK(data_transformed,model,rawID,'recon2');
%
for i=1:length(prot_sens)
    prot_sens{i,3} = Expression_Model(find(ismember(Expression_Model(:,1),prot_sens{i,2})),8);
end
% Now, prot_sens is a dataframe where data_transformed genes and associated
% model-genes are connected. Now. it is ready for splitting into EPI and
% MES. 
EPI_genes = num2cell(zeros(length(FDR_test),1));
EPI_genes = cat(2,EPI_genes,num2cell(num2cell(zeros(length(EPI_genes),1))));
EPI_genes(:,1) = num2cell(FDR_test(:));
MES_genes = EPI_genes;
for i=1:length(FDR_test)
    for k=1:length(prot_sens{i,3})
        if nanmean(cell2mat(Expression_Model(find(ismember(Expression_Model(:,8),prot_sens{i,3}(k))),2:4))) > nanmean(cell2mat(Expression_Model(find(ismember(Expression_Model(:,8),prot_sens{i,3}(k))),5:7)))
            EPI_genes{i,2}{k} = prot_sens{i,3}{k};
        elseif nanmean(cell2mat(Expression_Model(find(ismember(Expression_Model(:,8),prot_sens{i,3}(k))),2:4))) < nanmean(cell2mat(Expression_Model(find(ismember(Expression_Model(:,8),prot_sens{i,3}(k))),5:7)))
            MES_genes{i,2}{k} = prot_sens{i,3}{k};
        else
            continue
        end
    end
        EPI_genes{i,2} = EPI_genes{i,2}';
        MES_genes{i,2} = MES_genes{i,2}';
end

% Don't know why, but there is a 0 in the gene lists... Take it out!
for i=1:length(EPI_genes)
    for k=1:length(EPI_genes{i,2})
        if EPI_genes{i,2}{k} == 0
            EPI_genes{i,2}{k} = [];
        else
            EPI_genes{i,2}{k} = EPI_genes{i,2}{k};
        end
    end
end

for i=1:length(MES_genes)
    for k=1:length(MES_genes{i,2})
        if MES_genes{i,2}{k} == 0
            MES_genes{i,2}{k} = [];
        else
            MES_genes{i,2}{k} = MES_genes{i,2}{k};
        end
    end
end


for i=1:length(EPI_genes)
    EPI_genes{i,2} = EPI_genes{i,2}(~cellfun('isempty',EPI_genes{i,2}));
end

for i=1:length(MES_genes)
    MES_genes{i,2} = MES_genes{i,2}(~cellfun('isempty',MES_genes{i,2}));
end

% And now to account for genes that are only present in one cell line
% (excluded in t-test since NaN values are overlooked). 
EPI_count = zeros(length(Expression_Model),1);
MES_count = EPI_count;

for i=1:length(Expression_Model)
    EPI_count(i) = length(find(~isnan(cell2mat(Expression_Model(i,2:4)))));
    MES_count(i) = length(find(~isnan(cell2mat(Expression_Model(i,5:7)))));
end

EPI_and_MES_count = [num2cell(EPI_count) num2cell(MES_count) num2cell(zeros(length(EPI_count),1))];
idx = find(abs(cell2mat(EPI_and_MES_count(:,1)) - cell2mat(EPI_and_MES_count(:,2))) > 1);
idx2 = [num2cell(idx) num2cell(zeros(length(idx),1))];

for i=1:length(idx)
    if length(find(~isnan(cell2mat(Expression_Model(idx(i),2:4))))) > length(find(~isnan(cell2mat(Expression_Model(idx(i),5:7)))))
        idx2(i,2) = {'EPI'};
    else
        idx2(i,2) = {'MES'};
    end
end

EPI_only_genes = Expression_Model(idx(find(strcmp(idx2(:,2),'EPI')),1),8);
MES_only_genes = Expression_Model(idx(find(strcmp(idx2(:,2),'MES')),1),8);

for i=1:10
    EPI_genes{i,2} = cat(1,EPI_genes{i,2},EPI_only_genes);
    MES_genes{i,2} = cat(1,MES_genes{i,2},MES_only_genes);
end
% Use the FDR p-value < 0.05:
EPI_used = EPI_genes{5,2};
MES_used = MES_genes{5,2};

% Create the models:
tol = 1e-5;

% The amount of flux downregulation (1 corresponds to knockout)
downRegFrac=0.7;
downRegFrac_abs=0.7;

nrxns=length(model.rxns);
newlb=zeros(nrxns,1);
newub=zeros(nrxns,1);
lower=zeros(nrxns,1);
upper=zeros(nrxns,1);

for i=1:length(model.rxns)
    % Default values correspond to the original upper/lower bounds
    newlb(i)=model.lb(i);
    newub(i)=model.ub(i);
    if newlb(i) == newub(i)
        % Fixed reaction, possibly disabled
        fprintf('%d\t%s has fixed upper/lower bounds (%1.4f)\n', i, model.rxns{i}, model.lb(i))
        continue
    end
    flux=sampleMetaOut.points(i,:);
    lower(i)=min(flux);
    upper(i)=max(flux);
    newlb(i)=min(0,prctile(flux(flux<tol), downRegFrac*100));
    newub(i)=max(0,prctile(flux(flux>tol), (1 - downRegFrac)*100));
end

% Create the EPI and MES models
modelNew = model
modelNew.lb=lower;
modelNew.ub=upper;

% Constrain using microarray data:
idx_M=affectedRxns(model,EPI_used);
idx_E=affectedRxns(model,MES_used);

fprintf('EPI model: Number of affected reactions: %d\n', length(idx_E))
modelE=modelNew;
for i=1:length(idx_E)
    rxnE=idx_E(i);
    modelE=changeRxnBounds(modelE, modelE.rxns{rxnE}, newlb(rxnE),'l');
    modelE=changeRxnBounds(modelE, modelE.rxns{rxnE}, newub(rxnE),'u');
    fprintf('EPI\t%s\t%1.4f\t%1.4f\n', modelE.rxns{rxnE}, modelE.lb(rxnE), modelE.ub(rxnE))
end

for i=1:length(modelE.rxns)
    % Default values correspond to the original upper/lower bounds
    flux_E=sampleMetaOut.points(i,:);
    newlb_E(i)=min(0,prctile(flux_E(flux_E<tol), downRegFrac*100));
    newub_E(i)=max(0,prctile(flux_E(flux_E>tol), (1 - downRegFrac)*100));
end


fprintf('*************************************************\n')
fprintf('\nMES model: Number of affected reactions: %d\n', length(idx_M))
modelM=modelNew;
for i=1:length(idx_M)
    rxnM=idx_M(i);
    modelM=changeRxnBounds(modelM, modelM.rxns{rxnM}, newlb(rxnM),'l');
    modelM=changeRxnBounds(modelM, modelM.rxns{rxnM}, newub(rxnM),'u');
    fprintf('MES\t%s\t%1.4f\t%1.4f\n', modelM.rxns{rxnM}, modelM.lb(rxnM), modelM.ub(rxnM))
end

for i=1:length(modelM.rxns)
    % Default values correspond to the original upper/lower bounds
    flux_M = sampleMetaOut.points(i,:);
    newlb_M(i)=min(0,prctile(flux_M(flux_M<tol), downRegFrac*100));
    newub_M(i)=max(0,prctile(flux_M(flux_M>tol), (1 - downRegFrac)*100));
end

% Set all bounds to MS data
forcedExchangesE = [];
forcedExchangesM = [];
percentError = 0.1;
Ex = [];

for i = 2:size(fluxdataEpi0to48,1)
    Ex = MetNames2ExRxns((find(ismember(MetNames2ExRxns(:,1),fluxdataEpi0to48(i,1)))),2);
    if ~isempty(find(ismember(modelE.rxns,Ex)))
        
        forcedExchangesE = cat(1,forcedExchangesE,[Ex,cell2mat(fluxdataEpi0to48(i,2))-(abs(cell2mat(fluxdataEpi0to48(i,2))*percentError)),cell2mat(fluxdataEpi0to48(i,2))+(abs(cell2mat(fluxdataEpi0to48(i,2))*percentError))]);
        modelE=changeRxnBounds(modelE, Ex, cell2mat(fluxdataEpi0to48(i,2))-abs((cell2mat(fluxdataEpi0to48(i,2))*percentError)), 'l');
        modelE=changeRxnBounds(modelE, Ex, cell2mat(fluxdataEpi0to48(i,2))+abs((cell2mat(fluxdataEpi0to48(i,2))*percentError)), 'u');
       
        forcedExchangesM = cat(1,forcedExchangesM,[Ex,cell2mat(fluxdataMes0to48(i,2))-abs((cell2mat(fluxdataMes0to48(i,2))*percentError)),cell2mat(fluxdataMes0to48(i,2))+abs((cell2mat(fluxdataMes0to48(i,2))*percentError))]);
        modelM=changeRxnBounds(modelM, Ex, cell2mat(fluxdataMes0to48(i,2))-abs((cell2mat(fluxdataMes0to48(i,2))*percentError)), 'l');
        modelM=changeRxnBounds(modelM, Ex, cell2mat(fluxdataMes0to48(i,2))+abs((cell2mat(fluxdataMes0to48(i,2))*percentError)), 'u');

    else
    end
end

% To check for the metabolites which have huge errors in MS data
% Glutamine measured in a different experiment for accuracy
modelE.lb(find(ismember(modelE.rxns,'EX_pyr(e)'))) = -10;
modelE.ub(find(ismember(modelE.rxns,'EX_pyr(e)'))) = 0;
modelE.lb(find(ismember(modelE.rxns,'EX_lac_D(e)'))) = 0;
modelE.ub(find(ismember(modelE.rxns,'EX_lac_D(e)'))) = 100;

modelM.lb(find(ismember(modelM.rxns,'EX_pyr(e)'))) = -10;
modelM.ub(find(ismember(modelM.rxns,'EX_pyr(e)'))) = 0;
modelM.lb(find(ismember(modelM.rxns,'EX_lac_D(e)'))) = 0;
modelM.ub(find(ismember(modelM.rxns,'EX_lac_D(e)'))) = 100;
idx_relax =find(ismember(model.rxns,model.rxns));

% Relax models to achieve feasibility:
FBAE = optimizeCbModel(modelE);
FBAE.f
FBAM = optimizeCbModel(modelM);
FBAM.f
if FBAE.f == 0 || FBAM.f == 0
    try
    [modelENew,r,d,f,rxns_relaxed] = relax_rxns(modelE,idx_relax,0.01,[],[],0);
    [modelMNew,r,d,f,rxns_relaxed] = relax_rxns(modelM,idx_relax,0.01,[],[],0);
    catch
    end
else
    modelENew = modelE
    modelMNew = modelM
end

% ----------------- FBA -------------------------
% Check for feasibility of relaxed models:
FBAE = optimizeCbModel(modelENew);
if FBAE.f < tol
    fprintf('*************************************************\n')
    fprintf('Infeasible EPI model!\n')
else
    fprintf('*************************************************\n')
    fprintf('Feasible EPI model achieved!\n')
end

FBAM = optimizeCbModel(modelMNew);
if FBAM.f < tol
    fprintf('Infeasible MES model!\n')
else
    fprintf('Feasible MES model achieved!\n')
end

% ----------------- GEA -------------------------
[grRatioE,~,~,~,hasEffectE] = singleGeneDeletion(modelENew,'FBA',modelENew.genes,'true');
lethalE = find(isnan(grRatioE));
[grRatioM,~,~,~,hasEffectM] = singleGeneDeletion(modelMNew,'FBA',modelMNew.genes,'true');
lethalM = find(isnan(grRatioM));
lethal_onlyM = model.genes(setdiff(lethalM,lethalE));
lethal_onlyE = model.genes(setdiff(lethalE,lethalM));

save('proteomic_EMT_GSMMs.mat')
clear
%% 4) Creation of the RNAseq-constrained iBreast2886 GSMMs
% Load necessary measurements
% From Halldorsson et al.:
load('EMTmodel.mat');
load('MetabolomicsFluxEpi0to48.mat');
load('MetabolomicsFluxMes0to48.mat');
load('Met2EXrxns.mat');

% Other:
load('ID_conversion_gene_names_to_ENTREZ.mat');
load('sampled_iBreast2886.mat');

[~,~,raw_rna] = xlsread('D492_D492M_RNAseq.xlsx', 2); % Data from 
% Halldorsson et al., 2017. 
raw_rna2 = raw_rna(2:end,:);
raw_rna2(find(cellfun(@(x)ischar(x),raw_rna2(:,3))),:) = []; % Remove strings

[RNA_array_recon2] = ExpData_to_RxnNames_STK(raw_rna2,model,IDconv,'recon2');
genes_up = RNA_array_recon2(find(cell2mat(RNA_array_recon2(:,3)) > 0),5);
genes_down = RNA_array_recon2(find(cell2mat(RNA_array_recon2(:,3)) < 0),5);

% Create the models:
tol = 1e-5;

% The amount of flux downregulation (1 corresponds to knockout)
downRegFrac=0.7;
downRegFrac_abs=0.7;

nrxns=length(model.rxns);
newlb=zeros(nrxns,1);
newub=zeros(nrxns,1);
lower=zeros(nrxns,1);
upper=zeros(nrxns,1);

for i=1:length(model.rxns)
    % Default values correspond to the original upper/lower bounds
    newlb(i)=model.lb(i);
    newub(i)=model.ub(i);
    if newlb(i) == newub(i)
        % Fixed reaction, possibly disabled
        fprintf('%d\t%s has fixed upper/lower bounds (%1.4f)\n', i, model.rxns{i}, model.lb(i))
        continue
    end
    flux=sampleMetaOut.points(i,:);
    lower(i)=min(flux);
    upper(i)=max(flux);
    newlb(i)=min(0,prctile(flux(flux<tol), downRegFrac*100));
    newub(i)=max(0,prctile(flux(flux>tol), (1 - downRegFrac)*100));
end

% Create the EPI and MES models
modelNew = model
modelNew.lb=lower;
modelNew.ub=upper;

% Constrain using microarray data:
idx_M=affectedRxns(model,genes_down);
idx_E=affectedRxns(model,genes_up);

fprintf('EPI model: Number of affected reactions: %d\n', length(idx_E))
modelE=modelNew;
for i=1:length(idx_E)
    rxnE=idx_E(i);
    modelE=changeRxnBounds(modelE, modelE.rxns{rxnE}, newlb(rxnE),'l');
    modelE=changeRxnBounds(modelE, modelE.rxns{rxnE}, newub(rxnE),'u');
    fprintf('EPI\t%s\t%1.4f\t%1.4f\n', modelE.rxns{rxnE}, modelE.lb(rxnE), modelE.ub(rxnE))
end

for i=1:length(modelE.rxns)
    % Default values correspond to the original upper/lower bounds
    flux_E=sampleMetaOut.points(i,:);
    newlb_E(i)=min(0,prctile(flux_E(flux_E<tol), downRegFrac*100));
    newub_E(i)=max(0,prctile(flux_E(flux_E>tol), (1 - downRegFrac)*100));
end


fprintf('*************************************************\n')
fprintf('\nMES model: Number of affected reactions: %d\n', length(idx_M))
modelM=modelNew;
for i=1:length(idx_M)
    rxnM=idx_M(i);
    modelM=changeRxnBounds(modelM, modelM.rxns{rxnM}, newlb(rxnM),'l');
    modelM=changeRxnBounds(modelM, modelM.rxns{rxnM}, newub(rxnM),'u');
    fprintf('MES\t%s\t%1.4f\t%1.4f\n', modelM.rxns{rxnM}, modelM.lb(rxnM), modelM.ub(rxnM))
end

for i=1:length(modelM.rxns)
    % Default values correspond to the original upper/lower bounds
    flux_M = sampleMetaOut.points(i,:);
    newlb_M(i)=min(0,prctile(flux_M(flux_M<tol), downRegFrac*100));
    newub_M(i)=max(0,prctile(flux_M(flux_M>tol), (1 - downRegFrac)*100));
end

% Set all bounds to MS data
forcedExchangesE = [];
forcedExchangesM = [];
percentError = 0.1;
Ex = [];

for i = 2:size(fluxdataEpi0to48,1)
    Ex = MetNames2ExRxns((find(ismember(MetNames2ExRxns(:,1),fluxdataEpi0to48(i,1)))),2);
    if ~isempty(find(ismember(modelE.rxns,Ex)))
        
        forcedExchangesE = cat(1,forcedExchangesE,[Ex,cell2mat(fluxdataEpi0to48(i,2))-(abs(cell2mat(fluxdataEpi0to48(i,2))*percentError)),cell2mat(fluxdataEpi0to48(i,2))+(abs(cell2mat(fluxdataEpi0to48(i,2))*percentError))]);
        modelE=changeRxnBounds(modelE, Ex, cell2mat(fluxdataEpi0to48(i,2))-abs((cell2mat(fluxdataEpi0to48(i,2))*percentError)), 'l');
        modelE=changeRxnBounds(modelE, Ex, cell2mat(fluxdataEpi0to48(i,2))+abs((cell2mat(fluxdataEpi0to48(i,2))*percentError)), 'u');
       
        forcedExchangesM = cat(1,forcedExchangesM,[Ex,cell2mat(fluxdataMes0to48(i,2))-abs((cell2mat(fluxdataMes0to48(i,2))*percentError)),cell2mat(fluxdataMes0to48(i,2))+abs((cell2mat(fluxdataMes0to48(i,2))*percentError))]);
        modelM=changeRxnBounds(modelM, Ex, cell2mat(fluxdataMes0to48(i,2))-abs((cell2mat(fluxdataMes0to48(i,2))*percentError)), 'l');
        modelM=changeRxnBounds(modelM, Ex, cell2mat(fluxdataMes0to48(i,2))+abs((cell2mat(fluxdataMes0to48(i,2))*percentError)), 'u');

    else
    end
end

% To check for the metabolites which have huge errors in MS data
% Glutamine measured in a different experiment for accuracy
modelE.lb(find(ismember(modelE.rxns,'EX_pyr(e)'))) = -10;
modelE.ub(find(ismember(modelE.rxns,'EX_pyr(e)'))) = 0;
modelE.lb(find(ismember(modelE.rxns,'EX_lac_D(e)'))) = 0;
modelE.ub(find(ismember(modelE.rxns,'EX_lac_D(e)'))) = 100;

modelM.lb(find(ismember(modelM.rxns,'EX_pyr(e)'))) = -10;
modelM.ub(find(ismember(modelM.rxns,'EX_pyr(e)'))) = 0;
modelM.lb(find(ismember(modelM.rxns,'EX_lac_D(e)'))) = 0;
modelM.ub(find(ismember(modelM.rxns,'EX_lac_D(e)'))) = 100;
idx_relax =find(ismember(model.rxns,model.rxns));

% Relax models to achieve feasibility:
FBAE = optimizeCbModel(modelE);
FBAE.f
FBAM = optimizeCbModel(modelM);
FBAM.f
if FBAE.f == 0 || FBAM.f == 0
    try
    [modelENew,r,d,f,rxns_relaxed] = relax_rxns(modelE,idx_relax,0.01,[],[],0);
    [modelMNew,r,d,f,rxns_relaxed] = relax_rxns(modelM,idx_relax,0.01,[],[],0);
    catch
    end
else
    modelENew = modelE
    modelMNew = modelM
end

% ----------------- FBA -------------------------
% Check for feasibility of relaxed models:
FBAE = optimizeCbModel(modelENew);
if FBAE.f < tol
    fprintf('*************************************************\n')
    fprintf('Infeasible EPI model!\n')
else
    fprintf('*************************************************\n')
    fprintf('Feasible EPI model achieved!\n')
end

FBAM = optimizeCbModel(modelMNew);
if FBAM.f < tol
    fprintf('Infeasible MES model!\n')
else
    fprintf('Feasible MES model achieved!\n')
end

% ----------------- GEA -------------------------
[grRatioE,~,~,~,hasEffectE] = singleGeneDeletion(modelENew,'FBA',modelENew.genes,'true');
lethalE = find(isnan(grRatioE));
[grRatioM,~,~,~,hasEffectM] = singleGeneDeletion(modelMNew,'FBA',modelMNew.genes,'true');
lethalM = find(isnan(grRatioM));
lethal_onlyM = model.genes(setdiff(lethalM,lethalE));
lethal_onlyE = model.genes(setdiff(lethalE,lethalM));

save('rnaseq_EMT_GSMMs.mat')
clear