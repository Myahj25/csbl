
%doing flux variablity analysis on the ecoli model
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

%fva 
[minFlux,maxFlux]=fluxVariability(model);

BothFlux = [minFlux,maxFlux];

%creating a heatmap to compare min and max flux
xvalues = {'min','max'};
yvalues = model.rxns;
figure(1)
h5 = heatmap(xvalues,yvalues,BothFlux);

% just wanted to see what a regular plot would look like here
figure (2)
idk = plot(BothFlux);

