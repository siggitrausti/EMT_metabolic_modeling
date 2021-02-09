% Script to estimate the carbon contribution from glutamine into the
% reductive carboxylation pathway based on the randomly sampled data 
% from all the GSMMs

% Sigurdur Karvelsson

%% Prepare a matrix for results

% Random numbers are all sampled flux distributions within the sampling
% matrix (in this case 1 to 5800). 
random_numbers = [1:5800]';

% Also, define the transport reactions within the iBreast2886 model to add
% into the function. Decreases runtime of function substantially:
load('EMTmodel.mat');
transport_reactions = findTransRxns(model);

% Define a matrix for inputting the results from the PathwayActivity
% calculations systematically:
type_pat = cat(2,repelem({'Microarray'},10),repelem({'Proteomic'},10),repelem({'RNAseq'},10),repelem({'Media'},10));
cell_pat = repmat([{'EPI'},{'MES'}],1,20);
path_pat = repmat(repelem([{'Glutamine to citrate'},{'Glutamine to malate'},{'Glutamine to aspartate'},{'Glutamine to proline'},{'Glutamine to GSH'}],2),1,4);
nreplicates = 5800; % No. sampled flux distributions within each GSMM
mat = cat(1,type_pat,cell_pat,path_pat,num2cell(zeros(nreplicates,40)));

%% Calculate the PathwayActivity for all GSMMs
for i = 1:4
    models = i;
    k = (models-1)*10;
    PATHWAY_ACTIVITY_WORKING_SCRIPT
    mat(4:size(mat,1),k+1) = num2cell(percentageE1);
    mat(4:size(mat,1),k+2)= num2cell(percentageM1);
    mat(4:size(mat,1),k+3) = num2cell(percentageE2);
    mat(4:size(mat,1),k+4) = num2cell(percentageM2);
    mat(4:size(mat,1),k+5) = num2cell(percentageE3);
    mat(4:size(mat,1),k+6)= num2cell(percentageM3);
    mat(4:size(mat,1),k+7) = num2cell(percentageE4);
    mat(4:size(mat,1),k+8) = num2cell(percentageM4);
    mat(4:size(mat,1),k+9) = num2cell(percentageE5);
    mat(4:size(mat,1),k+10)= num2cell(percentageM5);
    clearvars -except mat random_numbers transport_reactions
end

%% Write to a .txt file
T = cell2table(mat)
writetable(T,'PathwayActivity_measure.txt','Delimiter','\t')
