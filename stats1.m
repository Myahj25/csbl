

%% aligning the flux vectors
% Get reaction lists from models
all_reactions = unique([CAF1model.rxns; CAF2model.rxns; NF1model.rxns; NF2model.rxns]);


CAF1flux_aligned = NaN(length(all_reactions), 1);
CAF2flux_aligned = NaN(length(all_reactions), 1);
NF1flux_aligned = NaN(length(all_reactions), 1);
NF2flux_aligned = NaN(length(all_reactions), 1);

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

%% Normalize flux vectors
growthRxnIndexC1 = findRxnIDs(CAF1model, 'MAR13082');
growthRxnIndexC2 = findRxnIDs(CAF2model, 'MAR13082'); 
growthRxnIndexN1 = findRxnIDs(NF1model, 'MAR13082'); 
growthRxnIndexN2 = findRxnIDs(NF2model, 'MAR13082'); 

CAF1flux_norm = CAF1flux_aligned / CAF1modelanalysis.v(growthRxnIndexC1);
CAF2flux_norm = CAF2flux_aligned / CAF2modelanalysis.v(growthRxnIndexC2);
NF1flux_norm = NF1flux_aligned / NF1modelanalysis.v(growthRxnIndexN1);
NF2flux_norm = NF2flux_aligned / NF2modelanalysis.v(growthRxnIndexN2);

%% Wilcoxon ranksum test
p_values = NaN(length(all_reactions), 1);  % Initialize p-values vector with NaN

% perform wilcoxon ranksum test for each reaction, skipping NaNs
for i = 1:length(all_reactions)
    % extract fluxes for the current reaction, ignoring NaNs
    NF_fluxes = [NF1flux_norm(i), NF2flux_norm(i)];
    CAF_fluxes = [CAF1flux_norm(i), CAF2flux_norm(i)];
    
    % check if there are any NaNs in the flux values
    if all(~isnan(NF_fluxes)) && all(~isnan(CAF_fluxes))
        % perform wilcoxon ranksum test only if both groups have no NaNs
        p_values(i) = ranksum(NF_fluxes, CAF_fluxes);
    else
        % skip this reaction if any NaN is found
        p_values(i) = NaN;
    end
end

%making a table of pvalues
pRxnTable = table(all_reactions,p_values);
