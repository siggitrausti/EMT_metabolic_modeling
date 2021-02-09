function [gene_names] = EntrezToGeneNames(entrez_genes,ID_conversion)
% Finally a function for fast conversion of Entrez IDs to gene names 
% Requires a gene-conversion array (ID_conversion_gene_names_to_ENTREZ.mat)

% Sigurdur Karvelsson, February 2019
%%
% 1. Take out the splice variant: 
for i=1:length(entrez_genes)
    dot = strfind(entrez_genes{i},'.');
    entrez_genes{i} = entrez_genes{i}(1:dot-1);
end
% 2. find the id of genes in ID_conversion array (gene names in column 1,
% entrez ids in column 2)
if isa(ID_conversion{1,2},'char') ~=1
    for i=1:length(ID_conversion)
    ID_conversion{i,2} = num2str(ID_conversion{i,2});
    end
end
gene_id = find(ismember(ID_conversion(:,2),entrez_genes));
gene_names = ID_conversion(gene_id,1);
gene_names = unique(gene_names);
end