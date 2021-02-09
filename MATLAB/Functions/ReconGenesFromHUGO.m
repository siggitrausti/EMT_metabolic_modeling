function [isogenes_all] = ReconGenesFromHUGO(model,IDconv,genes,ReconType)
%ReconGenesFromHUGO Summary of this function goes here
%   Basicly, this function is supposed to generate a list of model genes
%   (with suffix) from a list of HUGO genes (regular gene names)

% DEC 2020: Changed to include suffixes for Recon3D.



% Sigurdur Karvelsson, JUNE 2019


genes_found = IDconv(find(ismember(IDconv(:,1),genes)),2);
isogenes_all = [];
for k=1:length(genes_found)
    gene_of_interest = genes_found(k);
    gene_of_interest2 = cellstr(num2str(cell2mat(gene_of_interest)));
    %display(gene_of_interest);
    isogenes = [];
    for i=1:10
        nr = i;
        if strcmp(ReconType,'recon3d')
            str_add = strcat('_AT',num2str(nr));
        elseif strcmp(ReconType,'recon2')
            str_add = strcat('.',num2str(nr));
        end
        %display(cellstr(num2str(cell2mat(gene_of_interest))));
        %display(str_add);
        jojo = strcat(gene_of_interest2,str_add);
        %display(jojo);
        hello = find(contains(model.genes,jojo));
        if length(hello) ~=0
            isogenes = cat(1,isogenes,jojo);
            continue
        else 
            %display('No more isogenes found')
            break
        end
    end
    isogenes_all = cat(1,isogenes_all,isogenes);
end

end

