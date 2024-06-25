%% gene matlab exercise
%not completed still in progress
%1. loading the E.coli core model
fileName = 'ecoli_core_model.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end
%backward compatibility with primer requires relaxation of upper bound on ATPM
modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%setting the maximum glucose uptake rate to 18.5mmol/gDW/h
model = changeRxnBounds(model,'EX_glc(e)',-18.5,'l');

% allow for unlimited oxygen (aerobic conditions)  
model = changeRxnBounds(model,'EX_o2(e)',-1000,'l');

%optimizing for growth
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

%2. finding genes that start with G and T
genes = delRxns2;

geneG = '';
for i = 1:length(genes)
    if startsWith(genes{i}, 'G')
        geneG = genes{i};
        break;
    end
end

geneT = '';
for i = 1:length(genes)
    if startsWith(genes{i}, 'T')
        geneT = genes{i};
        break;
    end
end

disp(['Gene starting with G: ', geneG]);
disp(['Gene starting with T: ', geneT]);

%3. knocking out these reactions

model = changeRxnBounds(model,'GLNabc',0,'l');

model = changeRxnBounds(model,'THD2',0,'l');
growthRate = optimizeCbModel(model);

%4. try and obtain 10 different solutions that will reach the same growth rate


