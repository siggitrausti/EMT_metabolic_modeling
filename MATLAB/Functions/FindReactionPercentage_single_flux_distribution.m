function [rel_flux,rxn_used] = FindReactionPercentage_single_flux_distribution(model,met1,met2,flux_vector,transport_reactions)
% A function for identifying the total percentage of flux between two
% metabolites, i.e. the percentage of total flux coming from metabolite 1
% to metabolite 2, even though it is through multiple reactions...

% As opposed to FindReactionPercentage, this function actually just does a
% single flux distribution. The values from multiple flux distributions can
% then be merged to get some sort of a distribution of values...


% Sigurdur Karvelsson, DECEMBER 2020

    if ~strcmp(findMetaboliteLocation(model,met1,1),findMetaboliteLocation(model,met2,1))
        is_between_compartments = 1;
    else
        is_between_compartments = 0;
    end
    
    if is_between_compartments == 1
        rel_flux = 1;
        %error_outp = 0; % NOTE CHANGE TO GET APPROPRIATE VALUES...
        rxn_used = {'TRANSPORT REACTION - TREATED AS 1'}; 
    else
        met_ID = findMetIDs(model,met1);
        
        
        
        comp_TMP = [];
        %error_outp_TMP2 = [];
        rxn_used_TMP2 = [];
        rel_flux_outp2 = [];
        met_str = cell2mat(model.mets(met_ID));
        met_str = met_str(1:end-2);
        met_all = find(startsWith(model.mets,met_str)); % Identifies all locations of this metabolite
        %transRxns = findTransRxns(model); % SUPER SLOW

        for o = 1:length(met_all)
         % Just so that no values are actually zero
            met_ID = findMetIDs(model,model.mets(met_all(o)));
            met_1_tmp = model.mets(met_all(o));
            comp_temp = findMetaboliteLocation(model,met_1_tmp,1); % Here it would be appropriate to account for transport reactions
            met2_tmp1 = cell2mat(model.mets(find(ismember(model.mets,met2))));
            met2_tmp2 = strcat(met2_tmp1(1:end-3),comp_temp);
            reactions = findRxnsFromMets(model,met_1_tmp);

            % Take out transrxns if there are any:
            if length(reactions) > 1
                reactions2 = setdiff(reactions,transport_reactions);
                if ~isempty(reactions2) & is_between_compartments == 0 % This should also apply for every reaction, not just when the reactions are more than 1...
                    reactions = reactions2;
                end
            end
            rxn_IDs = zeros(length(reactions),1);
            for i=1:length(reactions)
                rxn_IDs(i) = findRxnIDs(model,reactions(i));
            end

            %prod_cons_cell = cell(length(reactions),3);
            %prod_cons_cell(:,1) = reactions;
            %prod_cons_cell(:,3) = printRxnFormula(model,reactions,false);
            flux_vec = flux_vector(rxn_IDs); % For actual flux values
            flux_vec_2 = zeros(size(flux_vec)); % For directionality values
            flux_vec_3 = zeros(size(flux_vec)); % For percentage values

            for i=1:length(reactions)
                if sign(model.S(met_ID,rxn_IDs(i)))*sign(flux_vec(i)) == 1 
                    flux_vec_2(i) = 1; % For producing reaction
                elseif sign(model.S(met_ID,rxn_IDs(i)))*sign(flux_vec(i)) == -1
                    flux_vec_2(i) = -1; % For consuming reaction
                end
            end
            id_prod = find(flux_vec_2(:)==1);
            id_cons = find(flux_vec_2(:)==-1);
            total_producing_flux = sum(abs(flux_vec(id_prod)));
            total_consuming_flux = sum(abs(flux_vec(id_cons)));
            flux_vec_3(id_prod) = (abs(flux_vec(id_prod))/total_producing_flux); % producing percentage added
            flux_vec_3(id_cons) = -(abs(flux_vec(id_cons))/total_consuming_flux); % consuming percentage added
              
            temp_rxns = {};
            for j=1:length(reactions)
                temp_mets = findMetsFromRxns(model,reactions(j));
                for l = 1:length(temp_mets)
                    if strcmp(temp_mets(l),met2_tmp2) ~= 0 & flux_vec_2(j) ~=1
                        temp_rxns = [temp_rxns,reactions(j)]; % This currently only identifies whether or not the
                        % reactions are in the same reacion, but not whether it is
                        % a producing or consuming reaction....
                    end
                end
            end
            temp_rxns = temp_rxns';

            % 
            id_rxns = find(ismember(reactions,temp_rxns));
            rxn_used_TMP = reactions(id_rxns);
            total_val = sum(flux_vec_3(id_rxns));
            rel_flux_tmp = abs(total_val);
            compartment = findMetaboliteLocation(model,model.mets(met_all(o)),1);
            %display('Compartment')
            %display(compartment)
            %display(percentage_outp_TMP)
            comp_TMP = [comp_TMP,cellstr(compartment)];
            rxn_used_TMP2 = cat(1,rxn_used_TMP2,rxn_used_TMP);
            rel_flux_outp2 = [rel_flux_outp2,rel_flux_tmp];
        end
        relative_cont = [];
        for i = 1:length(comp_TMP)
            if strcmp(cell2mat(comp_TMP(i)),'[c]') == 1
                relative_cont = [relative_cont,0.54];
            elseif strcmp(cell2mat(comp_TMP(i)),'[m]') == 1
                relative_cont = [relative_cont,0.22];
                elseif strcmp(cell2mat(comp_TMP(i)),'[r]') == 1
                relative_cont = [relative_cont,0.12];
                elseif strcmp(cell2mat(comp_TMP(i)),'[n]') == 1
                relative_cont = [relative_cont,0.06];
                elseif strcmp(cell2mat(comp_TMP(i)),'[g]') == 1
                relative_cont = [relative_cont,0.03];
                elseif strcmp(cell2mat(comp_TMP(i)),'[x]') == 1
                relative_cont = [relative_cont,0.01];
                elseif strcmp(cell2mat(comp_TMP(i)),'[l]') == 1
                relative_cont = [relative_cont,0.01];
                elseif strcmp(cell2mat(comp_TMP(i)),'[e]') == 1 % Dont know how to handle, but fuck it
                relative_cont = [relative_cont,1];
            end
        end
        relative_cont(find(rel_flux_outp2 == 0)) = [];
        rel_flux_outp2(find(rel_flux_outp2 == 0)) = [];
        %display(rel_flux_outp2);
        relative_cont = MakeMDV(relative_cont);
        %display(relative_cont);
        rel_flux = rel_flux_outp2.*relative_cont;
        rel_flux = sum(rel_flux);
        %error_outp = error_outp_TMP2.*relative_cont; % NOTE CHANGE TO GET APPROPRIATE VALUES...
        %error_outp = sqrt(sum(error_outp.^2)); % NOTE CHANGE TO GET APPROPRIATE VALUES...
        rxn_used = unique(rxn_used_TMP2);
    end


