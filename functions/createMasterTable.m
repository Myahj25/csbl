function ordered_master_table = createMasterTable(all_reactions, source_table1, source_table2, source_table3, source_table4,reactionNames)
    % Create the master table with all reactions and empty columns for flux values
    num_reactions = length(all_reactions);
    
    master_table = table(all_reactions, reactionNames, NaN(num_reactions, 1), NaN(num_reactions, 1), NaN(num_reactions, 1), NaN(num_reactions, 1), ...
                         'VariableNames', {'Reaction', 'Reaction_Name', 'CAF1_Flux', 'NF1_Flux', 'CAF2_Flux', 'NF2_Flux'});
   

    % Populate master table with fluxes from each source table
    master_table = mapFluxToMaster(master_table, source_table1, 'CAF1_Flux', 'NF1_Flux');  % Populates CAF1 and NF1
    master_table = mapFluxToMaster(master_table, source_table2, 'CAF2_Flux', 'NF2_Flux');  % Populates CAF2 and NF2

    % Identify common reactions across tables and move them to the top
    common_reactions = intersect(intersect(intersect(source_table1.Reaction, source_table2.Reaction), ...
                                           source_table3.Reaction), ...
                                 source_table4.Reaction);
    isCommon = ismember(master_table.Reaction, common_reactions);
    common_part = master_table(isCommon, :);
    non_common_part = master_table(~isCommon, :);
    ordered_master_table = [common_part; non_common_part];
end


function master_table = mapFluxToMaster(master_table, source_table, caf_col, nf_col)
    % Ensure source table has the expected columns
    if all(ismember({'Reaction', 'Flux1', 'Flux2'}, source_table.Properties.VariableNames))
        [isInSource, idxMaster] = ismember(source_table.Reaction, master_table.Reaction);
        master_table{idxMaster(isInSource), caf_col} = source_table.Flux1(isInSource);
        master_table{idxMaster(isInSource), nf_col} = source_table.Flux2(isInSource);
    else
        error('Source table is missing required columns: Reaction, Flux1, Flux2.');
    end
end
