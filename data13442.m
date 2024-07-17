

%importing expression data 
CAF1_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442CAF-P6-9-26-20_abundance.csv');
CAF2_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442CAF-P5-9-6-20_abundance.csv');
NF1_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442NF-P5-RNA-KG-2-12-21_abundance.csv');
NF2_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442NF-P4-RNA-KG-2-12-21_abundance.csv');

%making context specific models with ftINIT
CAF1_data_struct.genes = CAF1_data{:, 1}; % gene names
CAF1_data_struct.tissues = CAF1_data.Properties.VariableNames(2);
CAF1_data_struct.levels = CAF1_data{:, 2}; % gene TPM values
CAF1_data_struct.threshold = 1;

%CAFmodel1 = ftINIT(prepData, CAF1_data_struct.tissues{1}, [], [], CAF1_data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);
% second model
CAF2_data_struct.genes = CAF2_data{:, 1}; % gene names
CAF2_data_struct.tissues = CAF2_data.Properties.VariableNames(2);
CAF2_data_struct.levels = CAF2_data{:, 2}; % gene TPM values
CAF2_data_struct.threshold = 1;

%CAFmodel2 = ftINIT(prepData, CAF2_data_struct.tissues{1}, [], [], CAF2_data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);

%third model 
NF1_data_struct.genes = NF1_data{:, 1}; % gene names
NF1_data_struct.tissues = NF1_data.Properties.VariableNames(2);
NF1_data_struct.levels = NF1_data{:, 2}; % gene TPM values
NF1_data_struct.threshold = 1;
%NFmodel1 = ftINIT(prepData, NF1_data_struct.tissues{1}, [], [], NF1_data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);

%fourth model 
NF2_data_struct.genes = NF2_data{:, 1}; % gene names
NF2_data_struct.tissues = NF2_data.Properties.VariableNames(2);
NF2_data_struct.levels = NF2_data{:, 2}; % gene TPM values
NF2_data_struct.threshold = 1;

%NFmodel2 = ftINIT(prepData, NF2_data_struct.tissues{1}, [], [], NF2_data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);



% CAFmodel1
% CAFmodel2
% NFmodel1
% NFmodel2

CAFmodel1 = changeObjective(CAFmodel1,'MAR13082');

CAFFBA1 = solveLP(CAFmodel1)
COBCAFFBA1 = optimizeCbModel(CAFmodel1)
NFmodel1 = changeObjective(NFmodel1,'MAR13082');

NFFBA1 = solveLP(NFmodel1);

% CAF1fluxdata = CAFFBA1.x;
% NF1fluxdata = NFFBA1.x;
% CAF1HM = heatmap(CAF1fluxdata);
% NF1HM = heatmap(NF1fluxdata);
% 
% 

models1 = {CAFmodel1,NFmodel1};
 compAttempt = compareMultipleModels(models1);

% rxn2Dmap = tsne(compAttempt.reactions.matrix', 'Distance', 'hamming', 'NumDimensions', 2, 'Perplexity', 2);
% 
% % plot and label the GEMs in tSNE space
% scatter(rxn2Dmap(:,1), rxn2Dmap(:,2));
% hold on
% text(rxn2Dmap(:,1), rxn2Dmap(:,2), compAttempt.modelIDs);