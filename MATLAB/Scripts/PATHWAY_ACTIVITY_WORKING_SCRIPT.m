% PathwayActivity script - Define metabolite lists and call the
% PathwayActivity function to determine the fractional contribution of a
% metabolite of interest to other metabolites (e.g. to compare to isotope
% tracer analyses).

% Sigurdur Karvelsson

%% Load correct sampled models:
if models == 1
    load('sampled_microarray_models.mat');
    model = samplesE_32hours;
elseif models == 2
    load('sampled_proteomic_models.mat');
    model = samplesE_32hours;
elseif models == 3
    load('sampled_rnaseq_models.mat');
    model = samplesE_32hours;
elseif models == 4
    load('sampled_media_models.mat');
    model = samplesE_32hours;
end

%% Call the PathwayActivity function for predetermined series of reactions:

% location does not matter, as the PathwayActivity function quantifies the 
% pathway in all compartments simultaneously

% Reductive carboxylation (Glutamine to citrate)
met_list1 = [{'gln_L[e]'},{'gln_L[c]'},{'glu_L[c]'},{'glu_L[m]'},{'akg[m]'},{'icit[m]'},{'cit[m]'},{'cit[c]'}]';
[percentageE1] = PathwayActivity(model,samplesE_32hours.points(:,random_numbers),met_list1,transport_reactions);
[percentageM1] = PathwayActivity(model,samplesM_32hours.points(:,random_numbers),met_list1,transport_reactions);

% Glutamine to malate
met_list2 = [{'gln_L[e]'},{'gln_L[c]'},{'glu_L[c]'},{'glu_L[m]'},{'akg[m]'},{'icit[m]'},{'cit[m]'},{'cit[c]'},{'oaa[c]'},{'mal_L[c]'}]';
[percentageE2] = PathwayActivity(model,samplesE_32hours.points(:,random_numbers),met_list2,transport_reactions);
[percentageM2] = PathwayActivity(model,samplesM_32hours.points(:,random_numbers),met_list2,transport_reactions);

% Glutamine to aspartate
met_list3 = [{'gln_L[e]'},{'gln_L[c]'},{'glu_L[c]'},{'glu_L[m]'},{'akg[m]'},{'icit[m]'},{'cit[m]'},{'cit[c]'},{'oaa[c]'},{'asp_L[c]'}]';
[percentageE3] = PathwayActivity(model,samplesE_32hours.points(:,random_numbers),met_list3,transport_reactions);
[percentageM3] = PathwayActivity(model,samplesM_32hours.points(:,random_numbers),met_list3,transport_reactions);

% Glutamine to proline
met_list4 = [{'gln_L[e]'},{'gln_L[c]'},{'glu_L[c]'},{'glu_L[m]'},{'glu5sa[m]'},{'1pyr5c[m]'},{'pro_L[m]'}]';
[percentageE4] = PathwayActivity(model,samplesE_32hours.points(:,random_numbers),met_list4,transport_reactions);
[percentageM4] = PathwayActivity(model,samplesM_32hours.points(:,random_numbers),met_list4,transport_reactions);

% Glutamine to glutathione
met_list5 = [{'gln_L[e]'},{'gln_L[c]'},{'glu_L[c]'},{'glucys[c]'},{'gthrd[c]'}]';
[percentageE5] = PathwayActivity(model,samplesE_32hours.points(:,random_numbers),met_list5,transport_reactions);
[percentageM5] = PathwayActivity(model,samplesM_32hours.points(:,random_numbers),met_list5,transport_reactions);




