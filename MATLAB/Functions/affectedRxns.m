function [idxAffectedRxns]=affectedRxns(model, geneList)
% Returns a list of the reactions affected by the genes in geneList
%
% Based on the deleteModelGenes function.

feasible=true;
idxAffectedRxns=[];
% Find gene indices in model
[isInModel,geneInd] = ismember(geneList,model.genes);
if (any(isInModel))
   % If there are any zero elements in geneInd remove them from the 
   % geneList and the geneInd because they correspond to genes that 
   % are not in the model.
   geneInd = geneInd( find( geneInd ) );

   % Find rxns associated with this gene
   rxnInd = find(any(model.rxnGeneMat(:,geneInd),2));
   if (~isempty(rxnInd))
       x = true(size(model.genes));
       x(geneInd) = false;
       constrainRxn = false(length(rxnInd),1);
       % Figure out if any of the reaction states is changed
       for j = 1:length(rxnInd)
           if (~eval(model.rules{rxnInd(j)}))
               constrainRxn(j) = true;
           end
       end
       idxAffectedRxns=rxnInd(constrainRxn);
   end
end
