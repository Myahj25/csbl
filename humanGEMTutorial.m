

load(['Human-GEM.mat']);
ihuman

%By default, the model objective (defined by the .c model field) is set to maximize flux
% through the generic human biomass reaction (MAR13082), and all exchange reactions are open.
ihuman.rxns(ihuman.c == 1)
%running FBA
sol = solveLP(ihuman)
solvercobra = optimizeCbModel (ihuman)

%% running ftINIT section

%ftINIT requires a preparation step, I commented it out after I ran it the first time
%prepData = prepHumanModelForftINIT(ihuman, false, '/Users/myah/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt', '/Users/myah/Human-GEM/model/reactions.tsv');
%save('prepData.mat', 'prepData')

gtex_data = readtable('/Users/myah/Documents/MATLAB/gtexSampForTutorialTPM.txt');
[~, n] = size(gtex_data);
numSamp = n-2; %the first two columns are the genes in ENSEMBL and gene symbols format

% take a look at the first few rows and columns of the table
gtex_data(1:5, 1:5)

data_struct.genes = gtex_data{:, 1}; % gene names
data_struct.tissues = gtex_data.Properties.VariableNames(3:n); % sample (tissue) names
data_struct.levels = gtex_data{:, 3:n}; % gene TPM values
data_struct.threshold = 1;

%running ftINIT without the second step 
model1 = ftINIT(prepData, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);

model1

%running ftINIT with the second step
model2 = ftINIT(prepData, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), false, true);

model2

%It is recommended to change the model id to a more descriptive name
model1.id = data_struct.tissues{1};

%running ftINIT for all samples 
for i = 1:numSamp
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
    models{i} = ftINIT(prepData, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);
 end
 
 save('models.mat', 'models')

%% getting data to export to R to make visuals 
baseModel = prepData.refModel;

% now build a matrix saying which reactions are on
compMat = false(length(baseModel.rxns), length(models));

for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,models{i}.rxns);
end

% run t-sne
rng(1);  %set random seed to make reproducible
proj_coords = tsne(double(compMat.'), 'Distance', 'hamming', 'NumDimensions', 2, 'Exaggeration', 6, 'Perplexity', 10);

% export to R
d = struct();
d.tsneX = proj_coords(:, 1);
d.tsneY = proj_coords(:, 2);
save('TSNE.mat', 'd');
