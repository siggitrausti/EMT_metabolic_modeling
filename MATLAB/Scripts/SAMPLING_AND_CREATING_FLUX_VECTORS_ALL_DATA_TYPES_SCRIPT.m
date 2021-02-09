% Script to generate the sampling matrices and median flux vectors from 
% the sampled data from all the GSMMs

% Sigurdur Karvelsson

%% Sample all models (Takes 8 whole days):
steps = 5800 % steps for each point. 
time = 32*60*60; % 32 hours in seconds
% microarray
load('microarray_EMT_GSMMs.mat');
modelENew.csense(1:size(modelENew.S,1),1) = 'E';
modelMNew.csense(1:size(modelMNew.S,1),1) = 'E';
[samplesE_32hours, MF_E] = gpSampler(modelENew,steps,[],time);
[samplesM_32hours, MF_M] = gpSampler(modelMNew,steps,[],time);
clearvars -except samplesE_32hours samplesM_32hours modelENew modelMNew MF_E MF_M
save('sampled_microarray_models.mat');

% proteomic
load('proteomic_EMT_GSMMs.mat');
modelENew.csense(1:size(modelENew.S,1),1) = 'E';
modelMNew.csense(1:size(modelMNew.S,1),1) = 'E';
[samplesE_32hours, MF_E] = gpSampler(modelENew,steps,[],time);
[samplesM_32hours, MF_M] = gpSampler(modelMNew,steps,[],time);
clearvars -except samplesE_32hours samplesM_32hours modelENew modelMNew MF_E MF_M
save('sampled_proteomic_models.mat');

% rnaseq
load('rnaseq_EMT_GSMMs.mat');
modelENew.csense(1:size(modelENew.S,1),1) = 'E';
modelMNew.csense(1:size(modelMNew.S,1),1) = 'E';
[samplesE_32hours, MF_E] = gpSampler(modelENew,steps,[],time);
[samplesM_32hours, MF_M] = gpSampler(modelMNew,steps,[],time);
clearvars -except samplesE_32hours samplesM_32hours modelENew modelMNew MF_E MF_M
save('sampled_rnaseq_models.mat');

% media
load('media_EMT_GSMMs.mat');
modelENew.csense(1:size(modelENew.S,1),1) = 'E';
modelMNew.csense(1:size(modelMNew.S,1),1) = 'E';
[samplesE_32hours, MF_E] = gpSampler(modelENew,steps,[],time);
[samplesM_32hours, MF_M] = gpSampler(modelMNew,steps,[],time);
clearvars -except samplesE_32hours samplesM_32hours modelENew modelMNew MF_E MF_M
save('sampled_media_models.mat');


%% Create flux image for Escher. 

% For microarray data:
load('sampled_microarray_models.mat');
flux_E = FluxVectorize(modelENew,samplesE_32hours.points);
flux_M = FluxVectorize(modelMNew,samplesM_32hours.points);
Escher('microarray_fluxes_DEC2020',modelENew,flux_M,flux_E)

% For proteomic data:
load('sampled_proteomic_models.mat');
flux_E = FluxVectorize(modelENew,samplesE_32hours.points);
flux_M = FluxVectorize(modelMNew,samplesM_32hours.points);
Escher('proteomic_fluxes_DEC2020',modelENew,flux_M,flux_E)

% For RNAseq data:
load('sampled_rnaseq_models.mat');
flux_E = FluxVectorize(modelENew,samplesE_32hours.points);
flux_M = FluxVectorize(modelMNew,samplesM_32hours.points);
Escher('rnaseq_fluxes_DEC2020',modelENew,flux_M,flux_E)

% For media data:
load('sampled_media_models.mat');
flux_E = FluxVectorize(modelENew,samplesE_32hours.points);
flux_M = FluxVectorize(modelMNew,samplesM_32hours.points);
Escher('media_fluxes_DEC2020',modelENew,flux_M,flux_E)


