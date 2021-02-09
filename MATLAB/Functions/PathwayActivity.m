function [predicted_labelling_percentage_vector] = PathwayActivity(model,sampled_results,met_list,transport_reactions)

% A function that works on individual flux distributions to identify 
% the cumulative consumption of a product metabolite of interest. 

predicted_labelling_percentage_vector = zeros(size(sampled_results,2),1);
for j = 1:size(sampled_results,2)
    percentage_total = zeros(length(met_list)-1,1);
    for i=1:length(met_list)-1
        metabolite1 = met_list(i);
        metabolite2 = met_list(i+1);
        rel_flux = FindReactionPercentage_single_flux_distribution(model,metabolite1,metabolite2,sampled_results(:,j),transport_reactions); % TEST MARCH 2020
        percentage_total(i) = rel_flux;
    end
    % https://mattstats.wordpress.com/2013/04/19/numerical-stability/ 
    outp_pre = sum(log(percentage_total));
    predicted_labelling_percentage_TMP = outp_pre(end);
    predicted_labelling_percentage_vector(j) = predicted_labelling_percentage_TMP;
end