% Script for generating the breast cancer patients GSMMs from proteomic
% data from Tang et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6276229/

% Sigurdur Karvelsson

%% 1) Data processing:
initCobraToolbox;
tang = readtable('tang_proteomic_data_JAN2021.txt','ReadVariableNames',false,'Delimiter',',');
tang = table2cell(tang);
genes = readtable('tang_proteomic_genes_JAN2021.txt','ReadVariableNames',false);
genes = table2cell(genes);
genes = cat(1,[{'patient'},{'estrogen'}]',genes(:,2));
comb_dat = cat(2,genes,tang');
comb_dat = comb_dat([1,4:size(comb_dat,1)],:);

% Now create the expressionmat:
load('EMTmodel.mat');
load('ID_gene_names_to_ENTREZ.mat');
[Expression_Model] = ExpData_to_RxnNames_STK(comb_dat,model,IDconv,'recon2');

% Create the patient-based constraint matrix:
patient_names2 = tang(:,1)';% Proteomic
patient_names = cat(2,cat(2,{'genes'},patient_names2),{'entrez_recon'});
tcga = cat(1,patient_names,Expression_Model);
tcga_fin = tcga(:,2:end);
tcga_fin = cat(2,tcga_fin(:,end),tcga_fin);
tcga_fin = tcga_fin(:,1:end-1);
threshold_vector = zeros(length(tcga_fin),1);

% Get the list of constrained reactions within each patient:
for i=2:length(tcga_fin)
    threshold_vector(i) = prctile(cell2mat(tcga_fin(i,2:end)),60);
end

outp = zeros(size(tcga_fin));
for i=2:size(tcga_fin,2)
    for j=2:size(tcga_fin,1)
        if cell2mat(tcga_fin(j,i)) < threshold_vector(j)
            outp(j,i) = 1;
        else
            outp(j,i) = 0;
        end
    end
end

patient_vector = [patient_names2' num2cell(zeros(length(patient_names2'),1))]; % Proteomic
%
for i=1:size(patient_vector,1)
    temp_idx = find(outp(:,i+1) == 1);
    patient_vector{i,2} = tcga_fin(temp_idx,1);
end

%% 2) Create the patient-specific models:
load('sampled_iBreast2886.mat'); % Sampled iBreast2886 model for 32 hours
tol = 1e-5;
downRegFrac=0.7; % Same as with proteomic and microarray models
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
    %flux = samples_32hours.points(i,:);
    lower(i)=min(flux);
    upper(i)=max(flux);
    newlb(i)=min(0,prctile(flux(flux<tol), downRegFrac*100));
    newub(i)=max(0,prctile(flux(flux>tol), (1 - downRegFrac)*100));
end

% Create the EPI and MES models
modelNew = model;
modelNew.lb=lower;
modelNew.ub=upper;

% Pre-allocate vectors for RNAseq analysis:
lethal_rnaseq = num2cell(zeros(length(patient_vector),1));
FBA_relax_bool = num2cell(zeros(length(patient_vector),1));
FBA_flux = num2cell(zeros(length(patient_vector),1));
lethal_mets_ids = num2cell(zeros(length(patient_vector),1));
for q=1:length(patient_vector)
    idx_E=affectedRxns(model,patient_vector{q,2}); % have to put up an iterative algorithm... 

    fprintf('Model: %s\n',cell2mat(patient_vector(q)))
    fprintf('Number of affected reactions: %d\n', length(idx_E))
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
    idx_relax = find(ismember(model.rxns,model.rxns));
    % Relax models to achieve feasibility:
    FBAE = optimizeCbModel(modelE);
    FBAE.f

    if FBAE.f == 0
        try
        [modelENew,r,d,f,rxns_relaxed] = relax_rxns(modelE,idx_relax,0.001,[],[],0);
        catch 
            fprintf('Inconsistent data in iteration %d \n', q);
        end
        relax_bool = 1;
    else
        modelENew = modelE;
        relax_bool = 0;
    end

    if relax_bool == 1 & ~exist('modelENew')
        FBA_relax_bool{q} = relax_bool;
        FBA_flux{q} = FBAE.x;
        lethal_rnaseq{q} = 0;
        continue
    elseif relax_bool ==1 & exist('modelENew')
            % Check for feasibility of relaxed models:
        FBAE = optimizeCbModel(modelENew);
        if FBAE.f < tol
            fprintf('*************************************************\n')
            fprintf('Infeasible EPI model!\n')
        else
            fprintf('*************************************************\n')
            fprintf('Feasible EPI model achieved!\n')
        end

        [grRatioE,~,~,~,hasEffectE] = singleGeneDeletion(modelENew,'FBA',modelENew.genes,'true');
        lethalE = find(isnan(grRatioE));
        lethal_rnaseq{q} = lethalE;
        FBA_relax_bool{q} = relax_bool;
        FBA_flux{q} = FBAE.x;
        clear modelENew
    elseif relax_bool ==0 & exist('modelENew')
                % Check for feasibility of relaxed models:
        FBAE = optimizeCbModel(modelENew);
        if FBAE.f < tol
            fprintf('*************************************************\n')
            fprintf('Infeasible EPI model!\n')
        else
            fprintf('*************************************************\n')
            fprintf('Feasible EPI model achieved!\n')
        end

        [grRatioE,~,~,~,hasEffectE] = singleGeneDeletion(modelENew,'FBA',modelENew.genes,'false');
        lethalE = find(isnan(grRatioE));
        lethal_rnaseq{q} = lethalE;
        FBA_relax_bool{q} = relax_bool;
        FBA_flux{q} = FBAE.x;
        clear modelENew
    end
end

save('GEA_of_all_patients_proteomic_TANG_PAPER_with_FBA.mat');
%% 3) Analyse the data:
% Read in data to split the patients into tumor vs nontumor:
metadat = readtable('tang_proteomic_sample_info_JAN2021.txt','ReadVariableNames',false);
metadat = table2cell(metadat);

tumor_pats = metadat(find(contains(metadat(:,3),'Tumor')),:);
patient_vector = patient_vector(find(contains(metadat(:,3),'Tumor')),:);
lethal_rnaseq = lethal_rnaseq(find(contains(metadat(:,3),'Tumor')),:);

% Find genes that are sign. enriched as being essential for each subtype:
erneg_pats = tumor_pats(find(contains(tumor_pats(:,4),'ER-')),2);
erpos_pats = tumor_pats(find(contains(tumor_pats(:,4),'ER+')),2);

id_pos = find(ismember(patient_vector(:,1),erpos_pats));
id_neg = find(ismember(patient_vector(:,1),erneg_pats));

erneg = lethal_rnaseq(id_neg); 
erpos = lethal_rnaseq(id_pos);

lethal_rnaseq2 = lethal_rnaseq(union(id_neg,id_pos));
gene_essentiality_count = num2cell(zeros(length(model.genes),1));
for i = 1:length(model.genes)
    count_temp = 0;
    for j=1:length(lethal_rnaseq2)
        if (find(ismember(lethal_rnaseq2{j},i)))
            count_temp = count_temp + 1;
        end
    end
    gene_essentiality_count{i} = count_temp;
end

% Go through the subsets:   
erneg_essentiality_count = num2cell(zeros(length(model.genes),1));
for i = 1:length(model.genes)
    count_temp = 0;
    for j=1:length(erneg)
        if (find(ismember(erneg{j},i)))
            count_temp = count_temp + 1;
        end
    end
    erneg_essentiality_count{i} = count_temp;
end
       
erpos_essentiality_count = num2cell(zeros(length(model.genes),1));
for i = 1:length(model.genes)
    count_temp = 0;
    for j=1:length(erpos)
        if (find(ismember(erpos{j},i)))
            count_temp = count_temp + 1;
        end
    end
    erpos_essentiality_count{i} = count_temp;
end

% Perform empirical p-value calcuations. First, sample similar sized
% patient groups (n = 32 and 33) 1000 times to estimate the null distribution:

npats = 32; 
perm_val = 1000;
random_patient_pos = cell(length(model.genes),perm_val);
for k = 1:perm_val
    rand_id = randsample(length(patient_vector),npats,true);
    rand_pats = lethal_rnaseq(rand_id);
    rand_essentiality_count = num2cell(zeros(length(model.genes),1));
    for i = 1:length(model.genes)
        count_temp = 0;
        for j=1:length(rand_pats)
            if (find(ismember(rand_pats{j},i)))
                count_temp = count_temp + 1;
            end
        end
        rand_essentiality_count{i} = count_temp;
    end
    random_patient_pos(:,k) = rand_essentiality_count;
end   

npats = 33;
perm_val = 1000;
random_patient_neg = cell(length(model.genes),perm_val);
for k = 1:perm_val
    rand_id = randsample(length(patient_vector),npats,true);
    rand_pats = lethal_rnaseq(rand_id);
    rand_essentiality_count = num2cell(zeros(length(model.genes),1));
    for i = 1:length(model.genes)
        count_temp = 0;
        for j=1:length(rand_pats)
            if (find(ismember(rand_pats{j},i)))
                count_temp = count_temp + 1;
            end
        end
        rand_essentiality_count{i} = count_temp;
    end
    random_patient_neg(:,k) = rand_essentiality_count;
end 

% For each gene, calculate the empirical cumulative distribution function
% to get a p-value for the ER-neg and positive patients:
p_val_neg = zeros(length(erneg_essentiality_count),1);
for i = 1:length(p_val_neg)
    erneg_count_temp = cell2mat(erneg_essentiality_count(i));
    amount_equal_or_more = length(find(cell2mat(random_patient_neg(i,:)) >= erneg_count_temp));
    p_val_neg(i) = (amount_equal_or_more+1)/(perm_val+1);
end

p_val_pos = zeros(length(erpos_essentiality_count),1);
for i = 1:length(p_val_pos)
    erpos_count_temp = cell2mat(erpos_essentiality_count(i));
    amount_equal_or_more = length(find(cell2mat(random_patient_pos(i,:)) >= erpos_count_temp));
    p_val_pos(i) = (amount_equal_or_more+1)/(perm_val+1);
end

% Output genes to HGNC genes:
erneg_genes = EntrezToGeneNames(model.genes(find(p_val_neg < 0.05)),IDconv);
erpos_genes = EntrezToGeneNames(model.genes(find(p_val_pos < 0.05)),IDconv);
