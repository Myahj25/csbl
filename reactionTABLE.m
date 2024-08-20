
allReactions = humanModel.rxns;

reactionTable = table(allReactions, 'VariableNames', {'Reaction'});
reactionTable.Reaction = allReactions; 
reactionTable.CAFmodel1 = zeros(length(allReactions), 1);
reactionTable.CAFmodel2 = zeros(length(allReactions), 1);
reactionTable.NFmodel1 = zeros(length(allReactions), 1);
reactionTable.NFmodel2 = zeros(length(allReactions), 1);

for i = 1:length(allReactions)
    reaction = allReactions{i};
  
    reactionTable.CAFmodel1(i) = ismember(reaction, CAFmodel1.rxns);
    reactionTable.CAFmodel2(i) = ismember(reaction, CAFmodel2.rxns);
    reactionTable.NFmodel1(i) = ismember(reaction, NFmodel1.rxns);
    reactionTable.NFmodel2(i) = ismember(reaction, NFmodel2.rxns);
end

rxnFormulas = printRxnFormula(humanModel, 'rxnAbbrList', humanModel.rxns, 'printFlag', false, 'lineChangeFlag', true, 'metNameFlag', true);


reactionTable.Formula = rxnFormulas;
reactionTable.Name = humanModel.rxnNames;

reactionTable = movevars(reactionTable, 'Formula', 'Before', 'CAFmodel1');
reactionTable = movevars(reactionTable, 'Name', 'Before', 'Formula');


