%% load data to create models
load('exprData.mat');
%creating CAF1 data table
CAF1_data.genes = exprData.genes(:,1);
CAF1_data.tissue = exprData.tissues(1,1);
CAF1_data.levels = exprData.genes_tpm_scaled(:,1);
CAF1_data.threshold = exprData.threshold_75(1,1);

genesTable = table(CAF1_data.genes, 'VariableNames', {'Genes'});
tissueTable = table(repmat(CAF1_data.tissue, height(genesTable), 1), 'VariableNames', {'Tissue'});
levelsTable = table(CAF1_data.levels, 'VariableNames', {'Levels'});
thresholdTable = table(repmat(CAF1_data.threshold, height(genesTable), 1), 'VariableNames', {'Threshold'});


CAF1_datatable = [genesTable, tissueTable, levelsTable, thresholdTable];

%creating CAF2 data table
CAF2_data.genes = exprData.genes(:,1);
CAF2_data.tissue = exprData.tissues(2,1);
CAF2_data.levels = exprData.genes_tpm_scaled(:,2);
CAF2_data.threshold = exprData.threshold_75(1,2);

genesTable = table(CAF2_data.genes, 'VariableNames', {'Genes'});
tissueTable = table(repmat(CAF2_data.tissue, height(genesTable), 1), 'VariableNames', {'Tissue'});
levelsTable = table(CAF2_data.levels, 'VariableNames', {'Levels'});
thresholdTable = table(repmat(CAF2_data.threshold, height(genesTable), 1), 'VariableNames', {'Threshold'});


CAF2_datatable = [genesTable, tissueTable, levelsTable, thresholdTable];


% creating NF1 data table
NF1_data.genes = exprData.genes(:,1);
NF1_data.tissue = exprData.tissues(11,1);
NF1_data.levels = exprData.genes_tpm_scaled(:,11);
NF1_data.threshold = exprData.threshold_75(1,11);

genesTable = table(NF1_data.genes, 'VariableNames', {'Genes'});
tissueTable = table(repmat(NF1_data.tissue, height(genesTable), 1), 'VariableNames', {'Tissue'});
levelsTable = table(NF1_data.levels, 'VariableNames', {'Levels'});
thresholdTable = table(repmat(NF1_data.threshold, height(genesTable), 1), 'VariableNames', {'Threshold'});

NF1_datatable = [genesTable, tissueTable, levelsTable, thresholdTable];

% creating NF2 data Table

NF2_data.genes = exprData.genes(:,1);
NF2_data.tissue = exprData.tissues(12,1);
NF2_data.levels = exprData.genes_tpm_scaled(:,12);
NF2_data.threshold = exprData.threshold_75(1,12);

genesTable = table(NF2_data.genes, 'VariableNames', {'Genes'});
tissueTable = table(repmat(NF2_data.tissue, height(genesTable), 1), 'VariableNames', {'Tissue'});
levelsTable = table(NF2_data.levels, 'VariableNames', {'Levels'});
thresholdTable = table(repmat(NF2_data.threshold, height(genesTable), 1), 'VariableNames', {'Threshold'});

NF2_datatable = [genesTable, tissueTable, levelsTable, thresholdTable];

%% create models
%uses CreateMODEL function
CAF1model = CreateMODEL(CAF1_datatable,prepData);
CAF2model = CreateMODEL(CAF2_datatable,prepData);
NF1model = CreateMODEL(NF1_datatable,prepData);
NF2model = CreateMODEL(NF2_datatable,prepData);

%% analyze specific media models using pFBA
%uses analyzeModel function
% CAF1model
[CAF1modelanalysis,CAF1mediamodel] = analyzeModel(CAF1model, 'MAR13082');

% CAF2model
[CAF2modelanalysis,CAF2mediamodel] = analyzeModel(CAF2model, 'MAR13082');

%NF1model
[NF1modelanalysis,NF1mediamodel] = analyzeModel(NF1model, 'MAR13082');

%NF2model
[NF2modelanalysis,NF2mediamodel] = analyzeModel(NF2model, 'MAR13082');

tasks = parseTaskList('/Users/myah/Downloads/metabolicTasks_Essential.txt');

% Check if the model can perform the tasks
CAF1results = checkTasks(CAF1mediamodel,'/Users/myah/Downloads/metabolicTasks_Essential.txt');
CAF2results = checkTasks(CAF2mediamodel,'/Users/myah/Downloads/metabolicTasks_Essential.txt');
NF1results = checkTasks(NF1mediamodel,'/Users/myah/Downloads/metabolicTasks_Essential.txt');
NF2results = checkTasks(NF2mediamodel,'/Users/myah/Downloads/metabolicTasks_Essential.txt');


%% aligning the flux vectors

%saving flux vectors 
CAF1flux = CAF1modelanalysis.v(:,1);
CAF2flux = CAF2modelanalysis.v(:,1);
NF1flux = NF1modelanalysis.v(:,1);
NF2flux = NF2modelanalysis.v(:,1);

% Get reaction lists from models
all_reactions = unique([CAF1model.rxns; CAF2model.rxns; NF1model.rxns; NF2model.rxns]);


CAF1flux_aligned = zeros(length(all_reactions), 1);
CAF2flux_aligned = zeros(length(all_reactions), 1);
NF1flux_aligned = zeros(length(all_reactions), 1);
NF2flux_aligned = zeros(length(all_reactions), 1);

% fill in the aligned flux vectors based on common reactions
for i = 1:length(all_reactions)
    rxn = all_reactions{i};
    
    % For CAF1model
    rxnIndexCAF1 = findRxnIDs(CAF1model, rxn);
    if rxnIndexCAF1 > 0
        CAF1flux_aligned(i) = CAF1flux(rxnIndexCAF1);
    end
    
    % For CAF2model
    rxnIndexCAF2 = findRxnIDs(CAF2model, rxn);
    if rxnIndexCAF2 > 0
        CAF2flux_aligned(i) = CAF2flux(rxnIndexCAF2);
    end
    
    % For NF1model
    rxnIndexNF1 = findRxnIDs(NF1model, rxn);
    if rxnIndexNF1 > 0
        NF1flux_aligned(i) = NF1flux(rxnIndexNF1);
    end
    
    % For NF2model
    rxnIndexNF2 = findRxnIDs(NF2model, rxn);
    if rxnIndexNF2 > 0
        NF2flux_aligned(i) = NF2flux(rxnIndexNF2);
    end
end

%% putting data in the right format
%getting reaction names
reactionIDs = all_reactions; % reaction IDs of aligned reactions 

% Initialize cell array for storing reaction names
reactionNames = cell(size(reactionIDs));

% Loop through each reaction ID and get the corresponding reaction name
for i = 1:length(reactionIDs)
    % Find the index of the reaction ID in the model
    reactionIndex = find(strcmp(ihuman.rxns, reactionIDs{i}));
    
    % If a match is found, store the corresponding reaction name
    if ~isempty(reactionIndex)
        reactionNames{i} = ihuman.rxnNames{reactionIndex};
    end
end

%% ranking flux differences
% uses rankFluxes function & createMasterTable

% normalizing fluxes
growthRxnIndexC1 = findRxnIDs(CAF1model, 'MAR13082');
growthRxnIndexC2 = findRxnIDs(CAF2model, 'MAR13082'); 
growthRxnIndexN1 = findRxnIDs(NF1model, 'MAR13082'); 
growthRxnIndexN2 = findRxnIDs(NF2model, 'MAR13082'); 

CAF1flux_norm = CAF1flux_aligned / CAF1modelanalysis.v(growthRxnIndexC1);
CAF2flux_norm = CAF2flux_aligned / CAF2modelanalysis.v(growthRxnIndexC2);
NF1flux_norm = NF1flux_aligned / NF1modelanalysis.v(growthRxnIndexN1);
NF2flux_norm = NF2flux_aligned / NF2modelanalysis.v(growthRxnIndexN2);


%CAF1 and NF1 comparison 
[normCAF1NF1_resultstable, normC1N1_top_diff] = rankFluxes(CAF1flux_norm,NF1flux_norm,all_reactions);

%CAF1 and NF2 comparison

[normCAF1NF2_resultstable, normC1N2_top_diff] = rankFluxes(CAF1flux_norm,NF2flux_norm,all_reactions);

%CAF2 and NF1 comparison

[normCAF2NF1_resultstable, normC2N1_top_diff] = rankFluxes(CAF2flux_norm,NF1flux_norm,all_reactions);

%CAF2 and NF2 comparison

[normCAF2NF2_resultstable, normC2N2_top_diff] = rankFluxes(CAF2flux_norm,NF2flux_norm,all_reactions);

% Extract the reactions from each table
normreactions1 = normC1N1_top_diff.Reaction;
normreactions2 = normC1N2_top_diff.Reaction;
normreactions3 = normC2N1_top_diff.Reaction;
normreactions4 = normC2N2_top_diff.Reaction;


normcommon_reactions_all = intersect(normreactions1, intersect(normreactions2, intersect(normreactions3, normreactions4)));

normordered_master_table = createMasterTable(all_reactions, normC1N1_top_diff, normC2N2_top_diff, normC2N1_top_diff, normC1N2_top_diff,reactionNames);

normflux_values = table(all_reactions, CAF1flux_norm, CAF2flux_norm, NF1flux_norm, NF2flux_norm, ...
                    'VariableNames', {'Reaction', 'CAF1_Flux', 'CAF2_Flux', 'NF1_Flux', 'NF2_Flux'});

% loops through each reaction in the master table
for i = 1:height(normordered_master_table)
    reaction = normordered_master_table.Reaction{i}; 
    % find the corresponding flux values in the normflux_values table
    idx = strcmp(normflux_values.Reaction, reaction);  
    
    if any(idx)  % Check if reaction exists in flux_values
        % Replace NaN values with corresponding flux values
        if isnan(normordered_master_table.CAF1_Flux(i))
            normordered_master_table.CAF1_Flux(i) = normflux_values.CAF1_Flux(idx);  % Replace NaN
        end
        if isnan(normordered_master_table.CAF2_Flux(i))
            normordered_master_table.CAF2_Flux(i) = normflux_values.CAF2_Flux(idx);  % Replace NaN
        end
        if isnan(normordered_master_table.NF1_Flux(i))
            normordered_master_table.NF1_Flux(i) = normflux_values.NF1_Flux(idx);  % Replace NaN
        end
        if isnan(normordered_master_table.NF2_Flux(i))
            normordered_master_table.NF2_Flux(i) = normflux_values.NF2_Flux(idx);  % Replace NaN
        end
    end
end

filename = 'normordered_master_table.xlsx';

writetable(normordered_master_table, filename, 'Sheet', 1, 'WriteRowNames', true);

rxnFormulas = printRxnFormula(ihuman, 'rxnAbbrList', all_reactions, 'printFlag', false, 'lineChangeFlag', true, 'metNameFlag', true);

% adding reaction formula to table
reactionFormulaMap = containers.Map(all_reactions, rxnFormulas);
normordered_master_table.Formula = cell(height(normordered_master_table), 1);
for i = 1:height(normordered_master_table)
    reactionID = normordered_master_table.Reaction{i};
    if isKey(reactionFormulaMap, reactionID)
        normordered_master_table.Formula{i} = reactionFormulaMap(reactionID);
    else
        normordered_master_table.Formula{i} = ''; % Leave empty if no formula is found
    end
end


normordered_master_table = movevars(normordered_master_table, 'Formula', 'Before', 'CAF1_Flux');

filename = 'normordered_master_table.xlsx';

writetable(normordered_master_table, filename, 'Sheet', 1, 'WriteRowNames', true);

%adding subsystems to table

subreactionIDs = CAF1model.rxns;       
subsystems = CAF1model.subSystems; 

% Create a map of reaction IDs to subsystems
reactionSubsystemMap = containers.Map(subreactionIDs, subsystems);


normordered_master_table.Subsystem = cell(height(normordered_master_table), 1);

for i = 1:height(normordered_master_table)
    reactionID = normordered_master_table.Reaction{i};
    if isKey(reactionSubsystemMap, reactionID)
        normordered_master_table.Subsystem{i} = reactionSubsystemMap(reactionID); 
    else
        normordered_master_table.Subsystem{i} = ''; 
    end
end


normordered_master_table = movevars(normordered_master_table, 'Subsystem', 'Before', 'CAF1_Flux');

filename = 'normordered_master_table.xlsx';

writetable(normordered_master_table, filename, 'Sheet', 1, 'WriteRowNames', true);


%% creating a table to display gene rules 

reactionToGeneMap = containers.Map(CAF1model.rxns, CAF1model.grRules);


specificReactions = normordered_master_table.Reaction; 


genesForReactions = cell(length(specificReactions), 1);

% Loop through each specific reaction and retrieve the genes
for i = 1:length(specificReactions)
    reactionID = specificReactions{i};
    if isKey(reactionToGeneMap, reactionID)
        genesForReactions{i} = reactionToGeneMap(reactionID);
    else
        genesForReactions{i} = 'No gene rule available'; 
    end
end


reactionGeneTable = table(specificReactions, genesForReactions, ...
    'VariableNames', {'Reaction', 'GeneRule'});


%% ranking Fold Changes
%uses calculateFoldChanges function

normC1N1dataTable = table(reactionIDs, reactionNames, CAF1flux_norm, NF1flux_norm, ...
    'VariableNames', {'Reactions', 'ReactionNames', 'CAF1flux_norm', 'NF1flux_norm'});
normC1N2dataTable = table(reactionIDs, reactionNames, CAF1flux_norm, NF2flux_norm, ...
    'VariableNames', {'Reactions', 'ReactionNames', 'CAF1flux_norm', 'NF2flux_norm'}); 
normC2N1dataTable = table(reactionIDs, reactionNames, CAF2flux_norm, NF1flux_norm, ...
    'VariableNames', {'Reactions', 'ReactionNames', 'CAF2flux_norm', 'NF1flux_norm'});
normC2N2dataTable = table(reactionIDs, reactionNames, CAF2flux_norm, NF2flux_norm, ...
    'VariableNames', {'Reactions', 'ReactionNames', 'CAF2flux_norm', 'NF2flux_norm'});

[normC1N1results,normC1N1_sigFlux] = calculateFluxFoldChanges(normC1N1dataTable, 'CAF1flux_norm', 'NF1flux_norm');
[normC1N2results,normC1N2_sigFlux] = calculateFluxFoldChanges(normC1N2dataTable,'CAF1flux_norm', 'NF2flux_norm');
[normC2N1results,normC2N1_sigFlux] = calculateFluxFoldChanges(normC2N1dataTable,'CAF2flux_norm', 'NF1flux_norm');
[normC2N2results,normC2N2_sigFlux] = calculateFluxFoldChanges(normC2N2dataTable,'CAF2flux_norm', 'NF2flux_norm');


normFoldreactions1 = normC1N1_sigFlux.Reactions;
normFoldreactions2 = normC1N2_sigFlux.Reactions;
normFoldreactions3 = normC2N1_sigFlux.Reactions;
normFoldreactions4 = normC2N2_sigFlux.Reactions;

% Find common reactions
normFoldcommon_reactions_12 = intersect(normFoldreactions1, normFoldreactions2);
normFoldcommon_reactions_13 = intersect(normFoldreactions1, normFoldreactions3);
normFoldcommon_reactions_14 = intersect(normFoldreactions1, normFoldreactions4);
normFoldcommon_reactions_23 = intersect(normFoldreactions2, normFoldreactions3);
normFoldcommon_reactions_24 = intersect(normFoldreactions2, normFoldreactions4);
normFoldcommon_reactions_34 = intersect(normFoldreactions3, normFoldreactions4);

normFoldcommon_reactions_all = intersect(normFoldreactions1, intersect(normFoldreactions2, intersect(normFoldreactions3, normFoldreactions4)));



% the reactions both comparisons have in common: comparing Fold vs Difference calculation

Bothcommon = intersect(normFoldcommon_reactions_all,normcommon_reactions_all);

%% creating fold change table

% Create maps for reaction formulas and subsystems
reactionFormulaMap = containers.Map(reactionIDs, rxnFormulas);  
reactionSubsystemMap = containers.Map(subreactionIDs, CAF1model.grRules); 


commonReactions = normFoldcommon_reactions_all; % common reactions across all the fold comparisons

% Initialize variables to store data
reactionIDs_common = cell(length(commonReactions), 1);
reactionNames_common = cell(length(commonReactions), 1);
formulas_common = cell(length(commonReactions), 1);
subsystems_common = cell(length(commonReactions), 1);
foldChange_C1N1 = zeros(length(commonReactions), 1);
foldChange_C1N2 = zeros(length(commonReactions), 1);
foldChange_C2N1 = zeros(length(commonReactions), 1);
foldChange_C2N2 = zeros(length(commonReactions), 1);

for i = 1:length(commonReactions)
    reactionID = commonReactions{i};
    
    % store reaction ID and name
    reactionIDs_common{i} = reactionID;
    reactionNames_common{i} = reactionNames{strcmp(reactionIDs, reactionID)};  

    % get reaction formula and subsystem
    if isKey(reactionFormulaMap, reactionID)
        formulas_common{i} = reactionFormulaMap(reactionID);
    else
        formulas_common{i} = 'N/A';
    end
    
    if isKey(reactionSubsystemMap, reactionID)
        subsystems_common{i} = reactionSubsystemMap(reactionID);
    else
        subsystems_common{i} = 'N/A';
    end
    
    % get fold changes from each comparison table
    foldChange_C1N1(i) = normC1N1_sigFlux.Log2FoldChange(strcmp(normC1N1_sigFlux.Reactions, reactionID));
    foldChange_C1N2(i) = normC1N2_sigFlux.Log2FoldChange(strcmp(normC1N2_sigFlux.Reactions, reactionID));
    foldChange_C2N1(i) = normC2N1_sigFlux.Log2FoldChange(strcmp(normC2N1_sigFlux.Reactions, reactionID));
    foldChange_C2N2(i) = normC2N2_sigFlux.Log2FoldChange(strcmp(normC2N2_sigFlux.Reactions, reactionID));
end

% Create final table
normFoldTable = table(reactionIDs_common, reactionNames_common, formulas_common, subsystems_common, ...
                   foldChange_C1N1, foldChange_C1N2, foldChange_C2N1, foldChange_C2N2, ...
                   'VariableNames', {'ReactionID', 'ReactionName', 'Formula', 'Subsystem', ...
                                     'FoldChange_C1N1', 'FoldChange_C1N2', 'FoldChange_C2N1', 'FoldChange_C2N2'});

% Display the table or write to Excel
disp(normFoldTable);
writetable(normFoldTable, 'Common_Reactions_with_FoldChanges.xlsx');


%% knocking out significant reactions from fold changes and differences

combinedReactions = union(normFoldcommon_reactions_all,normcommon_reactions_all); %the reactions from both the fold change analysis and difference calculations

CAF1model_k = CAF1mediamodel;
CAF2model_k = CAF2mediamodel;
NF1model_k = NF1mediamodel;
NF2model_k = NF2mediamodel;

% Running singleRxnDeletion for CAF and NF models
[allCAF1_grRatio, allCAF1_grRateKO, allCAF1_grRateWT, allCAF1_hasEffect, allCAF1_delRxn,allCAF1fluxSolution] = singleRxnDeletion(CAF1model_k,'FBA' ,CAF1model_k.rxns);

[allCAF2_grRatio, allCAF2_grRateKO, allCAF2_grRateWT, allCAF2_hasEffect, allCAF2_delRxn,allCAF2fluxSolution] = singleRxnDeletion(CAF2model_k, 'FBA', CAF2model_k.rxns);

[allNF1_grRatio, allNF1_grRateKO, allNF1_grRateWT, allNF1_hasEffect, allNF1_delRxn,allNF1fluxSolution] = singleRxnDeletion(NF1model_k, 'FBA', NF1model_k.rxns);

[allNF2_grRatio, allNF2_grRateKO, allNF2_grRateWT, allNF2_hasEffect, allNF2_delRxn,allNF2fluxSolution] = singleRxnDeletion(NF2model_k, 'FBA', NF2model_k.rxns);



% find indices for each model
[~, CAF1_indices] = ismember(CAF1model_k.rxns, ihuman.rxns);
[~, CAF2_indices] = ismember(CAF2model_k.rxns, ihuman.rxns);
[~, NF1_indices] = ismember(NF1model_k.rxns, ihuman.rxns);
[~, NF2_indices] = ismember(NF2model_k.rxns, ihuman.rxns);


mappedCAF1fluxSolution = zeros(length(ihuman.rxns), length(ihuman.rxns));
mappedCAF2fluxSolution = zeros(length(ihuman.rxns), length(ihuman.rxns));
mappedNF1fluxSolution = zeros(length(ihuman.rxns), length(ihuman.rxns));
mappedNF2fluxSolution = zeros(length(ihuman.rxns), length(ihuman.rxns));

% map each model's flux solution matrix to the Human-GEM reaction space
for i = 1:length(CAF1_indices)
    if CAF1_indices(i) > 0  
        for j = 1:length(CAF1_indices)
            if CAF1_indices(j) > 0  
                mappedCAF1fluxSolution(CAF1_indices(i), CAF1_indices(j)) = allCAF1fluxSolution(i, j);
            end
        end
    end
end

for i = 1:length(CAF2_indices)
    if CAF2_indices(i) > 0
        for j = 1:length(CAF2_indices)
            if CAF2_indices(j) > 0
                mappedCAF2fluxSolution(CAF2_indices(i), CAF2_indices(j)) = allCAF2fluxSolution(i, j);
            end
        end
    end
end

for i = 1:length(NF1_indices)
    if NF1_indices(i) > 0
        for j = 1:length(NF1_indices)
            if NF1_indices(j) > 0
                mappedNF1fluxSolution(NF1_indices(i), NF1_indices(j)) = allNF1fluxSolution(i, j);
            end
        end
    end
end

for i = 1:length(NF2_indices)
    if NF2_indices(i) > 0
        for j = 1:length(NF2_indices)
            if NF2_indices(j) > 0
                mappedNF2fluxSolution(NF2_indices(i), NF2_indices(j)) = allNF2fluxSolution(i, j);
            end
        end
    end
end

% saving
save('mappedCAF1fluxSolution.mat', 'mappedCAF1fluxSolution');
save('mappedCAF2fluxSolution.mat', 'mappedCAF2fluxSolution');
save('mappedNF1fluxSolution.mat', 'mappedNF1fluxSolution');
save('mappedNF2fluxSolution.mat', 'mappedNF2fluxSolution');


%% checking to see if i did it right
disp(['Number of non-zero entries in mappedCAF1fluxSolution: ', num2str(nnz(mappedCAF1fluxSolution))]);
disp(['Number of non-zero entries in allCAF1fluxSolution: ', num2str(nnz(allCAF1fluxSolution))]);


%% removing the reactions that have the same flux values


threshold = 0.1;  %% set a threshold of 0.1 to include reactions with small differences to get removed


sameFluxCAF1_NF1 = abs(mappedCAF1fluxSolution - mappedNF1fluxSolution) <= threshold;
sameFluxCAF1_NF2 = abs(mappedCAF1fluxSolution - mappedNF2fluxSolution) <= threshold;
sameFluxCAF2_NF1 = abs(mappedCAF2fluxSolution - mappedNF1fluxSolution) <= threshold;
sameFluxCAF2_NF2 = abs(mappedCAF2fluxSolution - mappedNF2fluxSolution) <= threshold;

% combine them
finalRemoveMask = sameFluxCAF1_NF1 | sameFluxCAF1_NF2 | sameFluxCAF2_NF1 | sameFluxCAF2_NF2;


% identify indices to keep
keepIndices = find(~any(finalRemoveMask, 1));


% filter the flux solutions
filteredCAF1fluxSolution = mappedCAF1fluxSolution(keepIndices, keepIndices);
filteredCAF2fluxSolution = mappedCAF2fluxSolution(keepIndices, keepIndices);
filteredNF1fluxSolution = mappedNF1fluxSolution(keepIndices, keepIndices);
filteredNF2fluxSolution = mappedNF2fluxSolution(keepIndices, keepIndices);

% check
disp(['Filtered CAF1 matrix size: ', num2str(size(filteredCAF1fluxSolution, 1)), ' x ', num2str(size(filteredCAF1fluxSolution, 2))]);
