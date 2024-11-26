function [results_table, top_diff_reactions] = rankFluxes(flux1_aligned, flux2_aligned, reactions)

    
    rank1 = tiedrank(flux1_aligned);  
    rank2 = tiedrank(flux2_aligned);  

    
    results_table = table(reactions, flux1_aligned, flux2_aligned, rank1, rank2, ...
                          'VariableNames', {'Reaction', 'Flux1', 'Flux2', 'Rank1', 'Rank2'});
    

    disp('Full Results Table:');
    disp(results_table);
    
    flux_difference = abs(flux1_aligned - flux2_aligned);


    [~, diff_rank] = sort(flux_difference, 'descend');

 
    results_table.Difference = flux_difference;  
    results_table.Diff_Rank = tiedrank(-flux_difference); 

  
   

   percentile_99 = prctile(flux_difference, 99);

    % Select reactions in the 99th percentile of flux differences
    top_diff_reactions = results_table(flux_difference >= percentile_99, :);

    % Sort the selected reactions by their flux difference rank
    top_diff_reactions = sortrows(top_diff_reactions, 'Diff_Rank');

    % Display top reactions in the 99th percentile for flux differences
    disp('Reactions in the 99th Percentile of Flux Differences:');
    disp(top_diff_reactions);
end
