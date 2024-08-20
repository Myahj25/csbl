

%importing expression data 
CAF1_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442CAF-P6-9-26-20_abundance.csv');
CAF2_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442CAF-P5-9-6-20_abundance.csv');
NF1_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442NF-P5-RNA-KG-2-12-21_abundance.csv');
NF2_data = readtable('/Users/myah/Documents/MATLAB/humangem_13442NF-P4-RNA-KG-2-12-21_abundance.csv');

% creating models
CAFmodel1 = CreateMODEL(CAF1_data,prepData);
CAFmodel2 = CreateMODEL(CAF2_data,prepData);
NFmodel1 = CreateMODEL(NF1_data,prepData);
NFmodel2 = CreateMODEL(NF2_data,prepData);

%comparing number of reaction
 [RxnNumCAF1,~] = size(CAFmodel1.rxns);
 [RxnNumCAF2,~] = size(CAFmodel2.rxns);
 [RxnNumNF1,~] = size(NFmodel1.rxns);
 [RxnNumNF2,~] = size(NFmodel2.rxns);


xALLrxns = categorical({'CAF model 1' ,'CAF model 2','NF model 1', 'NF model 2'});
yALLrxns = [RxnNumCAF1 RxnNumCAF2 RxnNumNF1 RxnNumNF2];
figure(1)
bar(xALLrxns,yALLrxns);
title('Number of Reactions in the Models','FontSize', 14);
ylabel('Reactions')

%comparing number of genes
[originalnumCAF1genes,~] = size(CAF1_data);
[newnumCAF1genes,~] = size(CAFmodel1.genes);
[newnumCAF2genes,~] = size(CAFmodel2.genes);
[newnumNF1genes,~] = size(NFmodel1.genes);
[newnumNF2genes,~] = size(NFmodel2.genes);

figure(2)
xallgenes = categorical({'Original Number of Genes' ,'CAF model 1','CAF model 2','NF model 1','NF model 2'});
yallgenes = [originalnumCAF1genes newnumCAF1genes newnumCAF2genes newnumNF1genes newnumNF2genes];
bar(xallgenes,yallgenes);
title('Orignal Number of Genes vs Number of Genes in Models','FontSize', 14);
ylabel('Genes');


%% common and unique reaction analysis

%CAF vs CAF
% common reactions
rxnCompCAF1CAF2 = intersect(CAFmodel1.rxns,CAFmodel2.rxns);

% different reactions
rxnDiffCAF1CAF2 = setxor(CAFmodel1.rxns,CAFmodel2.rxns);
CAF1uniqueCAF2RxnDiff = setdiff(CAFmodel1.rxns,CAFmodel2.rxns);
CAF2uniqueCAF1RxnDiff = setdiff(CAFmodel2.rxns,CAFmodel1.rxns);

CC1common = size(rxnCompCAF1CAF2);
CC1unique = size(CAF1uniqueCAF2RxnDiff);

CC2common = size(rxnCompCAF1CAF2);
CC2unique = size(CAF2uniqueCAF1RxnDiff);

xCCrxns = categorical({'CAF1 rxns' 'CAF2 rxns'});
yCCrxns = [CC1common CC1unique; CC2common CC2unique]; 
figure (3)
bar(xCCrxns,yCCrxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of CAF1 and CAF2','FontSize',14);


%NF vs NF
% common reactions
rxnCompNF1NF2 = intersect(NFmodel1.rxns,NFmodel2.rxns);

% different reactions
rxnDiffNF1NF2 = setxor(NFmodel1.rxns,NFmodel2.rxns);
NF1uniqueNF2RxnDiff = setdiff(NFmodel1.rxns,NFmodel2.rxns);
NF2uniqueNF1RxnDiff = setdiff(NFmodel2.rxns,NFmodel1.rxns);

NN1common = size(rxnCompNF1NF2);
NN1unique = size(NF1uniqueNF2RxnDiff);

NN2common = size(rxnCompNF1NF2);
NN2unique = size(NF2uniqueNF1RxnDiff);    

xNNrxns = categorical({'NF1 rxns' 'NF2 rxns'});
yNNrxns = [NN1common NN1unique; NN2common NN2unique]; 
figure (4)
bar(xNNrxns,yNNrxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of NF1 and NF2','FontSize',14);

%CAF1 vs NF1
% common reactions
rxnCompCAF1NF1 = intersect(CAFmodel1.rxns,NFmodel1.rxns);

% different reactions
rxnDiffCAF1NF1 = setxor(CAFmodel1.rxns,NFmodel1.rxns);
CAF1uniqueNF1RxnDiff = setdiff(CAFmodel1.rxns,NFmodel1.rxns);
NF1uniqueCAF1RxnDiff = setdiff(NFmodel1.rxns,CAFmodel1.rxns);

N1C1common = size(rxnCompCAF1NF1);
N1CAF1unique = size(CAF1uniqueNF1RxnDiff);

C1N1common = size(rxnCompCAF1NF1);
C1NF1unique = size(NF1uniqueCAF1RxnDiff);

xN1C1rxns = categorical({'CAF1 rxns' 'NF1 rxns'});
yN1C1rxns = [N1C1common N1CAF1unique; C1N1common C1NF1unique]; 
figure (5)
bar(xN1C1rxns,yN1C1rxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of NF1 and CAF1','FontSize',14);


%CAF1 vs NF2
% common reactions
rxnCompCAF1NF2 = intersect(CAFmodel1.rxns,NFmodel2.rxns);

% different reactions
rxnDiffCAF1NF2 = setxor(CAFmodel1.rxns,NFmodel2.rxns);
CAF1uniqueNF2RxnDiff = setdiff(CAFmodel1.rxns,NFmodel2.rxns);
NF2uniqueCAF1RxnDiff = setdiff(NFmodel2.rxns,CAFmodel1.rxns);

N2C1common = size(rxnCompCAF1NF2);
N2CAF1unique = size(CAF1uniqueNF2RxnDiff);

C1N2common = size(rxnCompCAF1NF2);
C1NF2unique = size(NF2uniqueCAF1RxnDiff);

xN2C1rxns = categorical({'CAF1 rxns' 'NF2 rxns'});
yN2C1rxns = [N2C1common N2CAF1unique; C1N2common C1NF2unique]; 
figure (6)
bar(xN2C1rxns,yN2C1rxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of NF2 and CAF1','FontSize',14);

%CAF2 vs NF1
% common reactions
rxnCompCAF2NF1 = intersect(CAFmodel2.rxns,NFmodel1.rxns);

% different reactions
rxnDiffCAF2NF1 = setxor(CAFmodel2.rxns,NFmodel1.rxns);
CAF2uniqueNF1RxnDiff = setdiff(CAFmodel2.rxns,NFmodel1.rxns);
NF1uniqueCAF2RxnDiff = setdiff(NFmodel1.rxns,CAFmodel2.rxns);

N1C2common = size(rxnCompCAF2NF1);
N1CAF2unique = size(CAF2uniqueNF1RxnDiff);

C2N1common = size(rxnCompCAF2NF1);
C2NF1unique = size(NF1uniqueCAF2RxnDiff);

xN1C2rxns = categorical({'CAF2 rxns' 'NF1 rxns'});
yN1C2rxns = [N1C2common N1CAF2unique; C2N1common C2NF1unique]; 
figure (7)
bar(xN1C2rxns,yN1C2rxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of NF1 and CAF2','FontSize',14);



%CAF2 vs NF2 
% common reactions
rxnCompCAF2NF2 = intersect(CAFmodel2.rxns,NFmodel2.rxns);

% different reactions
rxnDiffCAF2NF2 = setxor(CAFmodel2.rxns,NFmodel2.rxns);
CAF2uniqueNF2RxnDiff = setdiff(CAFmodel2.rxns,NFmodel2.rxns);
NF2uniqueCAF2RxnDiff = setdiff(NFmodel2.rxns,CAFmodel2.rxns);

N2C2common = size(rxnCompCAF2NF2);
N2CAF2unique = size(CAF2uniqueNF2RxnDiff);

C2N2common = size(rxnCompCAF2NF2);
C2NF2unique = size(NF2uniqueCAF2RxnDiff);

xN2C2rxns = categorical({'CAF2 rxns' 'NF2 rxns'});
yN2C2rxns = [N2C2common N2CAF2unique; C2N2common C2NF2unique]; 
figure (8)
bar(xN2C2rxns,yN2C2rxns,"stacked");
legend('common','unique');
title('Comparing Common and Unique Reactions of NF2 and CAF2','FontSize',14);

% prctile(CAF1_data{:, 2},75)
% prctile(CAF2_data{:, 2},75)
% 
% prctile(NF1_data{:, 2},75)
% prctile(NF2_data{:, 2},75)
% 
% temporaryData = [CAF1_data{:, 2}, CAF2_data{:, 2}, NF1_data{:, 2}, NF2_data{:, 2}];
% 
% figure(9)
% boxplot(temporaryData)

%rxnFormulas = printRxnFormula(humanModel, 'rxnAbbrList', humanModel.rxns(1:5), 'printFlag', 1, 'lineChangeFlag', 1, 'metNameFlag', 1);
