%importing expression data 
CAF1_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442CAF-P6-9-26-20_abundance.csv');
CAF2_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442CAF-P5-9-6-20_abundance.csv');
NF1_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442NF-P5-RNA-KG-2-12-21_abundance.csv');
NF2_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442NF-P4-RNA-KG-2-12-21_abundance.csv');

%% making context specific models with ftINIT
%first model 
CAF1_data_struct.genes = CAF1_data{:, 1}; % gene names
CAF1_data_struct.tissues = CAF1_data.Properties.VariableNames(2);
CAF1_data_struct.levels = CAF1_data{:, 2}; % gene TPM values
CAF1_data_struct.threshold = 26;

CAFmodel1 = ftINIT(prepData, CAF1_data_struct.tissues{1}, [], [], CAF1_data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);
% second model
CAF2_data_struct.genes = CAF2_data{:, 1}; % gene names
CAF2_data_struct.tissues = CAF2_data.Properties.VariableNames(2);
CAF2_data_struct.levels = CAF2_data{:, 2}; % gene TPM values
CAF2_data_struct.threshold = 27;

CAFmodel2 = ftINIT(prepData, CAF2_data_struct.tissues{1}, [], [], CAF2_data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);

%third model 
NF1_data_struct.genes = NF1_data{:, 1}; % gene names
NF1_data_struct.tissues = NF1_data.Properties.VariableNames(2);
NF1_data_struct.levels = NF1_data{:, 2}; % gene TPM values
NF1_data_struct.threshold = 29;
NFmodel1 = ftINIT(prepData, NF1_data_struct.tissues{1}, [], [], NF1_data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);

%fourth model 
NF2_data_struct.genes = NF2_data{:, 1}; % gene names
NF2_data_struct.tissues = NF2_data.Properties.VariableNames(2);
NF2_data_struct.levels = NF2_data{:, 2}; % gene TPM values
NF2_data_struct.threshold = 32;

NFmodel2 = ftINIT(prepData, NF2_data_struct.tissues{1}, [], [], NF2_data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);



% CAFmodel1
% CAFmodel2
% NFmodel1
% NFmodel2


%% FBA 
% CAFmodel1 = changeObjective(CAFmodel1,'MAR13082');
% 
% CAFFBA1 = solveLP(CAFmodel1)
% COBCAFFBA1 = optimizeCbModel(CAFmodel1)
% NFmodel1 = changeObjective(NFmodel1,'MAR13082');
% 
% NFFBA1 = solveLP(NFmodel1);

% CAF1fluxdata = CAFFBA1.x;
% NF1fluxdata = NFFBA1.x;


%% comparing models 
models1 = {CAFmodel1,NFmodel1};
 %compAttempt = compareMultipleModels(models1);

 %here I want to compare the number of reactions in each model 
 [RxnNumCAF1,~] = size(CAFmodel1.rxns);
 [RxnNumCAF2,~] = size(CAFmodel2.rxns);
 [RxnNumNF1,~] = size(NFmodel1.rxns);
 [RxnNumNF2,~] = size(NFmodel2.rxns);

xCAF = categorical({'Number of Rxns in CAF model 1' ,'Number of Rxns in CAF model 2'});
yCAF = [RxnNumCAF1, RxnNumCAF2];
figure(1)
bar(xCAF,yCAF);

xNF = categorical({'Number of Rxns in NF model 1' ,'Number of Rxns in NF model 2'});
yNF = [RxnNumNF1 RxnNumNF2];
figure(2)
bar(xNF,yNF)

xALL = categorical({'Number of Rxns in CAF model 1' ,'Number of Rxns in CAF model 2','Number of Rxns in NF model 1', 'Number of Rxns in NF model 2'});
yALL = [RxnNumCAF1 RxnNumCAF2 RxnNumNF1 RxnNumNF2];
figure(3)
bar(xALL,yALL);

%comparing number of genes
%CAF1
[originalnumCAF1genes,~] = size(CAF1_data);
[newnumCAF1genes,~] = size(CAFmodel1.genes);
figure(4)
xCAF1genes = categorical({'Original Number of Genes in CAF model 1' ,'New Number of Genes in CAF model 1'});
yCAF1genes = [originalnumCAF1genes newnumCAF1genes];
bar(xCAF1genes,yCAF1genes);

%CAF2
[originalnumCAF2genes,~] = size(CAF2_data);
[newnumCAF2genes,~] = size(CAFmodel2.genes);
figure(5)
xCAF2genes = categorical({'Original Number of Genes in CAF model 2' ,'New Number of Genes in CAF model 2'});
yCAF2genes = [originalnumCAF2genes newnumCAF2genes];
bar(xCAF2genes,yCAF2genes);

%NF1 
[originalnumNF1genes,~] = size(NF1_data);
[newnumNF1genes,~] = size(NFmodel1.genes);
figure(6)
xNF1genes = categorical({'Original Number of Genes in NF model 1' ,'New Number of Genes in NF model 1'});
yNF1genes = [originalnumNF1genes newnumNF1genes];
bar(xNF1genes,yNF1genes);

%NF2
[originalnumNF2genes,~] = size(NF2_data);
[newnumNF2genes,~] = size(NFmodel2.genes);
figure(7)
xNF2genes = categorical({'Original Number of Genes in NF model 2' ,'New Number of Genes in NF model 2'});
yNF2genes = [originalnumNF2genes newnumNF2genes];
bar(xNF2genes,yNF2genes);


%similar genes
geneCompCAF1 = intersect(CAFmodel1.genes,CAF1_data_struct.genes);
geneCompCAF2 = intersect(CAFmodel2.genes,CAF2_data_struct.genes);
geneCompNF1 = intersect(NFmodel1.genes,NF1_data_struct.genes);
geneCompNF2 = intersect(NFmodel2.genes,NF2_data_struct.genes);
geneCompCAF1NF1 = intersect(CAFmodel1.genes,NFmodel1.genes);
geneCompCAF2NF2 = intersect(CAFmodel2.genes,NFmodel2.genes);

%different genes
geneDiffCAF1 = setxor(CAFmodel1.genes,CAF1_data_struct.genes);
geneDiffCAF2 = setxor(CAFmodel2.genes,CAF2_data_struct.genes);
geneDiffNF1 = setxor(NFmodel1.genes,NF1_data_struct.genes);
geneDiffNF2 = setxor(NFmodel2.genes,NF2_data_struct.genes);
geneDiffCAF1NF1 = setxor(CAFmodel1.genes,NFmodel1.genes);
geneDiffCAF2NF2 = setxor(CAFmodel2.genes,NFmodel2.genes);

trialgeneDiffCAF1NF1 = setdiff(CAFmodel1.genes,NFmodel1.genes);

%similar rxns
rxnCompCAF1CAF2 = intersect(CAFmodel1.rxns,CAFmodel2.rxns);
rxnCompNF1NF2 = intersect(NFmodel1.rxns,NFmodel2.rxns);
rxnCompCAF1NF1 = intersect(CAFmodel1.rxns,NFmodel1.rxns);
rxnCompCAF2NF2 = intersect(CAFmodel2.rxns,NFmodel2.rxns);

%different rxn
rxnDiffCAF1CAF2 = setxor(CAFmodel1.rxns,CAFmodel2.rxns);
rxnDiffNF1NF2 = setxor(NFmodel1.rxns,NFmodel2.rxns);
rxnDiffCAF1NF1 = setxor(CAFmodel1.rxns,NFmodel1.rxns);
rxnDiffCAF2NF2 = setxor(CAFmodel2.rxns,NFmodel2.rxns);

trialrxnDiffCAF1NF1 = setdiff(CAFmodel1.rxns,NFmodel1.rxns);
trialrxnDiffNF1CAF1 = setdiff(NFmodel1.rxns,CAFmodel1.rxns);
trialrxnDiffCAF2NF2 = setdiff(CAFmodel2.rxns,NFmodel2.rxns);
trialrxnDiffNF2CAF2 = setdiff(NFmodel2.rxns,CAFmodel2.rxns);

%testing to make sure functions work
settrial = [1 3 4 5];
settrial2 = [1 3 5 6];
setres =setdiff(settrial2,settrial);
setres2 = setdiff(settrial,settrial2);
setres3=setxor(settrial,settrial2);

genetrial = '';
for i = 1:length(CAF1_data_struct.genes)
    if startsWith(CAF1_data_struct.genes{i}, 'ENSG00000087085')
        genetrial = CAF1_data_struct.genes{i};
        break;
    end
end

trialCAF1CAF2RxnDiff = setxor(trialrxnDiffCAF1NF1,trialrxnDiffCAF2NF2);

%% figures 

%common and unique factors 
%figure for CAF1 unique reactions within NF1 and vice versa 
CAF1unique = 133;
CAF1common = 8234;

NF1unique = 225;
NF1common = 8234; 

xCAF1NF1rxns = categorical({'CAF1 rxns' 'NF1 rxns'});
yCAF1NF1rxns = [CAF1common CAF1unique; NF1common NF1unique]; 
figure (8)
bar(xCAF1NF1rxns,yCAF1NF1rxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of CAF1 and NF1');

%figure for CAF2 unique reactions within NF2 and vice versa 
CAF2unique = 111;
CAF2common = 8220;

NF2unique = 250;
NF2common = 8220;

xCAF2NF2rxns = categorical({'CAF2 rxns' 'NF2 rxns'});
yCAF2NF2rxns = [CAF2common CAF2unique; NF2common NF2unique]; 
figure (9)
bar(xCAF2NF2rxns,yCAF2NF2rxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of CAF2 and NF2');


%% CAF1 and CAF2 comparisons
% common reactions
rxnCompCAF1CAF2 = intersect(CAFmodel1.rxns,CAFmodel2.rxns);

% different reactions
rxnDiffCAF1CAF2 = setxor(CAFmodel1.rxns,CAFmodel2.rxns);
CAF1CAF2RxnDiff = setdiff(CAFmodel1.rxns,CAFmodel2.rxns);
CAF2CAF1RxnDiff = setdiff(CAFmodel2.rxns,CAFmodel1.rxns);

%CAF1 common and unique reactions to CAF2 and vice versa
CC1common = 8240;
CC1unique = 127;

CC2common = 8240;
CC2unique = 91;

xCCrxns = categorical({'CAF1 rxns' 'CAF2 rxns'});
yCCrxns = [CC1common CC1unique; CC2common CC2unique]; 
figure (10)
bar(xCCrxns,yCCrxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of CAF1 and CAF2');


%% NF1 and NF2 comparisons
% common reactions
rxnCompNF1NF2 = intersect(NFmodel1.rxns,NFmodel2.rxns);

% different reactions
rxnDiffNF1NF2 = setxor(NFmodel1.rxns,NFmodel2.rxns);
NF1NF2RxnDiff = setdiff(NFmodel1.rxns,NFmodel2.rxns);
NF2NF1RxnDiff = setdiff(NFmodel2.rxns,NFmodel1.rxns);

%NF1 common and unique reactions to NF2 and vice versa 
NN1common = 8356;
NN1unique = 103;

NN2common = 8356;
NN2unique = 114;    

xNNrxns = categorical({'NF1 rxns' 'NF2 rxns'});
yNNrxns = [NN1common NN1unique; NN2common NN2unique]; 
figure (11)
bar(xNNrxns,yNNrxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of NF1 and NF2');

%comparing caf and nf (choosing CAF1 and NF2 because they have more unique
%reactions) 
% common reactions
rxnCompCAF1NF2 = intersect(CAFmodel1.rxns,NFmodel2.rxns);

% different reactions
rxnDiffCAF1NF2 = setxor(CAFmodel1.rxns,NFmodel2.rxns);
CAF1NF2RxnDiff = setdiff(CAFmodel1.rxns,NFmodel2.rxns);
NF2CAF1RxnDiff = setdiff(NFmodel2.rxns,CAFmodel1.rxns);

%NF1 common and unique reactions to NF2 and vice versa 
C1N2common = 8243;
C1N2unique = 124;

N2C1common = 8243;
N2C1unique = 227;    

xN2C1rxns = categorical({'CAF rxns' 'NF rxns'});
yN2C1rxns = [C1N2common C1N2unique; N2C1common N2C1unique]; 
figure (12)
bar(xN2C1rxns,yN2C1rxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of CAF and NF');

%comparing unique reactions with expressions profile 
CAF1unique_reactions = CAF1NF2RxnDiff;
CAF1all_reactions = CAFmodel1.rxns;
CAF1gene_data_table = CAFmodel1.grRules;
% Initialize a cell array to hold the indices
indices = cell(size(CAF1unique_reactions));

% Loop through each unique reaction
for i = 1:length(CAF1unique_reactions)
    % Find the indices of the current unique reaction in all_reactions
    indices{i} = find(strcmp(CAF1all_reactions, CAF1unique_reactions{i}));
end

% Display the gene data using indices
for i = 1:length(CAF1unique_reactions)
    disp(['Gene data for ', CAF1unique_reactions{i}, ':']);
    current_indices = indices{i}; % Get the indices for the current unique reaction
    gene_data = CAF1gene_data_table(current_indices); % Extract the gene data using the indices
    disp(gene_data); % Display the gene data
end
