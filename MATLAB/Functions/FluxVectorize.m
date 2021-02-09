function [flux_dist] = FluxVectorize(model,sampled_results,type)
% This function takes the results from random sampling, finds the median
% value of the histograms for each reaction of the model, and outputs a
% vector with a single value (the median) for every reaction. 
% Particularly useful for the 'MetabolicSummary' function and other
% FBA-based results (such as MOMA etc). 

if nargin < 3
    type = 'median'
end

flux_dist = zeros(length(model.rxns),1);
for i=1:length(flux_dist)
    if strcmp(type,'mean')
        flux_dist(i) = mean(sampled_results.points(i,:));
    elseif strcmp(type,'median')
        flux_dist(i) = median(sampled_results.points(i,:));
    end
end
