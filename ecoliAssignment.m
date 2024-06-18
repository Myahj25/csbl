
%% doing single reaction/gene deletion for aerobic conditions
%loading the E.coli core model
fileName = 'ecoli_core_model.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end
%backward compatibility with primer requires relaxation of upper bound on
%ATPM
modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%setting the maximum glucose uptake rate to 18.5mmol/gDW/h
model = changeRxnBounds(model,'EX_glc(e)',-18.5,'l');

% allow for unlimited oxygen 
model = changeRxnBounds(model,'EX_o2(e)',-1000,'l');

%setting biomass reaction as the objective function 
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

%single reaction deletion
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model);

xvalues = model.rxns;
yvalues = model.rxns;
figure(1)
h = heatmap(xvalues,yvalues,fluxSolution);
h.Title = 'Heatmap for single reaction deletion for aerobic conditions';

%single gene deletion
[grRatio2, grRateKO2, grRateWT2, hasEffect2, delRxns2, fluxSolution2] = singleGeneDeletion(model);


xvalues = model.genes;
yvalues = model.rxns;
figure(2)
h2 = heatmap(xvalues,yvalues,fluxSolution2);
h2.Title = 'Heatmap for single gene deletion for aerobic conditions';


%% doing single reaction/gene deletion for anaerobic conditions

model = changeRxnBounds(model,'EX_o2(e)',0,'l');

%setting biomass reaction as the objective function 
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

%single reaction deletion
[grRatio3, grRateKO3, grRateWT3, hasEffect3, delRxn3, fluxSolution3] = singleRxnDeletion(model);

xvalues = model.rxns; 
yvalues = model.rxns;
figure(3)
h3=heatmap(xvalues,yvalues,fluxSolution3);
h3.Title = 'Heatmap for single reaction deletion for anaerobic conditions';

%single gene deletion
[grRatio4, grRateKO4, grRateWT4, hasEffect4, delRxns4, fluxSolution4] = singleGeneDeletion(model);

xvalues = model.genes;
yvalues = model.rxns;
figure(4)
h4 = heatmap(xvalues,yvalues,fluxSolution4);
h4.Title = 'Heatmap for single gene deletion for anaerobic conditions';

