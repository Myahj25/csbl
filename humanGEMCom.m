
%comparing GEM models 
load('/Users/myah/Documents/MATLAB/tINIT_GEMs/run_tINIT_outputs/GTEx/tINIT_GTEx_outputs.mat');
INIT_output

model_ids = INIT_output.id;
models = INIT_output.model;

% this compares two or more condition-specific models generated from the same
% base model using high-dimensional comparisons in the reaction-space
res = compareMultipleModels(models);
res

% creating a clustergram to visualize the hamming distance
% GEMS that share fewer reactions willbe seperated by a larger Hamming
% distance and GEMS that contain many of the same reactions will have a
% smaller Hamming distance

clustergram(res.structComp, 'Symmetric', false, 'Colormap', 'bone', 'RowLabels', res.modelIDs, 'ColumnLabels', res.modelIDs);

% plot and label the GEMs in tSNE space
rxn2Dmap = tsne(res.reactions.matrix', 'Distance', 'hamming', 'NumDimensions', 2, 'Perplexity', 5);

scatter(rxn2Dmap(:,1), rxn2Dmap(:,2));
hold on
text(rxn2Dmap(:,1), rxn2Dmap(:,2), res.modelIDs);