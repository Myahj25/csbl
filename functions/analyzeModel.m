

function [modelAnalysis,outputmodel] = analyzeModel(inputModel, objectiveRxnID)
    
    % Keep a copy of the original model
    origModel = inputModel;

    % Identify and close uptake reactions
    [selExc, selUpt] = findExcRxns(inputModel);
    uptakeRxnIndices = find(selUpt);
    inputModel.lb(uptakeRxnIndices) = 0;

    % Change specific reaction bounds to -1000 for certain metabolites
    media = getMediaComponents("humangem");
    metabolites = media.humangem; 
    metsInModel = inputModel.metNames;  
    commonMets = intersect(metabolites, metsInModel); 

    metaboliteNames = metabolites;
    metaboliteIndices = [];
    metaboliteIDList = {};

    for i = 1:length(metaboliteNames)
        idx = find(strcmp(inputModel.metNames, metaboliteNames{i}));
        if ~isempty(idx)
            metaboliteIndices = [metaboliteIndices; idx];
            metaboliteIDList = [metaboliteIDList; inputModel.mets{idx}]; 
        end
    end

    metaboliteIndices = metaboliteIndices(inputModel.metComps(metaboliteIndices) == 1);

    
    metaboliteReactions = find(sum(abs(inputModel.S(metaboliteIndices,:)), 1) == 1)'; 

  
    [~, ia, ~] = intersect(metaboliteReactions, uptakeRxnIndices, 'legacy');
    rxnIds = metaboliteReactions(ia); 
    inputModel.lb(rxnIds) = -1000;  

   
    inputModel = changeObjective(inputModel, objectiveRxnID);
    modelAnalysis = optimizeCbModel(inputModel, 'max', 'one');
    outputmodel = inputModel;
end
