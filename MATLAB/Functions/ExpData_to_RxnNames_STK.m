function [fin_dat] = ExpData_to_RxnNames_STK(expData,model,IDconv,ReconType)

% An extension of previous ExpData_to_RxnNames by adding recon3d:
for i=1:size(expData,1)
    id_matches = ReconGenesFromHUGO(model,IDconv,expData(i,1),ReconType);
    %display(id_matches);
    if ~isempty(length(id_matches))
        outp_mat = repmat(cat(2,expData(i,:),num2cell(0)),length(id_matches),1);
        outp_mat(1:length(id_matches),end) = id_matches;
    else 
        continue
    end
    %display(outp_mat)
    if i == 1
        fin_dat = outp_mat;
    else
        fin_dat = cat(1,fin_dat,outp_mat);
    end
    %k = size(fin_dat,1);
    %display(fin_dat);
    
end