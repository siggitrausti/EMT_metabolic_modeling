% Write escher reaction files into filename
% Need model and absolute mean values of samplepoints.
% Or just any single flux data point for each reaction.

function Escher(filename,model, flux1, flux2)
if nargin < 3
    fprintf('You need more inputs stupid!\n');
    return
elseif nargin < 4 
    bSingleStage = true;
else
    bSingleStage = false;
end

% Some reaction names are different in Escher, need to address this!
for j = 1:length(flux1)
    if strcmpi(model.rxns{j}, 'EX_lac_L(e)') == 1
        model.rxns{j} = 'EX_lac__L_e'
    end
    if strcmpi(model.rxns{j}, 'GLCt2_2') == 1
        model.rxns{j} = 'GLCt1'
    end
     if strcmpi(model.rxns{j}, 'GNDc') == 1
        model.rxns{j} = 'GND'
    end
    if strcmpi(model.rxns{j}, 'EX_co2(e)') == 1
        model.rxns{j} = 'EX_co2_e'
    end
end

if bSingleStage
    if length(flux1) ~= length(model.rxns)
        fprintf('The number of fluxes and reactions must match!\n');
        return
    end
   
    
    fID = fopen(filename, 'wt'); 
    assert(fID ~= -1); % Assert stops the program if the specific condition is not
                       % met. Here fID = -1 means that the file could not be
                       % created.
    fprintf(fID, 'ID,stage\n');
    % now do a for loop that writes the reaction names and values for
    % every reaction.... let Escher worry about if we need them or not.
    for i = 1:length(flux1)
    fprintf(fID,'%s,%1.4f\n', model.rxns{i}, flux1(i)); 
    end
    fclose(fID);
    fprintf('Load %s into Escher for single-stage analysis :)\n', filename);
else % For the two stage analysis:
    
    % If this problem pops up it can possibly be resolved with
    % the FixEscher function.
    if length(flux1) ~= length(flux2)
        fprintf('The number of fluxes have to match!\n');
        return
    end
   
    if length(flux1) ~= length(model.rxns)
        fprintf('The number of fluxes and reactions must match!\n');
        return
    end
    
    fID = fopen(filename, 'wt'); 
    assert(fID ~= -1);
    fprintf(fID, 'ID,stage1,stage2\n');
    for i = 1:length(flux1)
    fprintf(fID,'%s,%1.4f,%1.4f\n', model.rxns{i}, flux1(i), flux2(i));
    end
    fclose(fID);
     fprintf('Load %s into Escher for two-stage analysis :)\n', filename);
end
