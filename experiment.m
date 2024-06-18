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

% allow for unlimited oxygen (this represents aerobic conditions
model = changeRxnBounds(model,'EX_o2(e)',-1000,'l');

%setting biomass reaction as the objective function 
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

%performing FBA with maximization of the biomass reaction as the objective
FBAsolution = optimizeCbModel(model,'max');

% this makes the flux distribution more convenient to analyze (if you want to
%see the flux distribution you can just uncomment out the print flux vector
%line below)
fluxData = FBAsolution.v;
nonZeroFlag = 1;
%printFluxVector(model, fluxData, nonZeroFlag)

%display an optimal flux vector on a metabolic map
outputFormatOK = changeCbMapOutput('matlab');
map=readCbMap('ecoli_core_map');
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;

%figure(1)
drawFlux(map,model,FBAsolution.v,options);

%simulating anaerobic conditions so change oxygen bounds 
model = changeRxnBounds(model,'EX_o2(e)',0,'l');

%calculating anaerobic growth rates
FBAsolution2 = optimizeCbModel(model,'max');

%comparing the flux vectors for both aerobic and anaerobic
fluxData2 = [FBAsolution.v,FBAsolution2.v];
%nonZeroFlag = 1;
%excFlag = 1;
printFluxVector(model, fluxData2);

%map for anaerobic conditions (in order to see this map you have to comment
%out figure 1)
figure(2)
%drawFlux(map, model, FBAsolution2.v, options);

%examining succinate
%this was a part of the tutorial but not the main part that i foucsed on so i commented 
% it out for clear results

% model = modelOri;
% model = changeRxnBounds(model,'EX_glc(e)',0,'l');
% model = changeRxnBounds(model,'EX_succ(e)',-20,'l');
% checkObjective(model);
% FBAsolution3 = optimizeCbModel(model,'max');
% FBAsolution3.f
% model = changeRxnBounds(model,'EX_o2(e)',0,'l');
% FBAsolution4 = optimizeCbModel(model,'max');
% FBAsolution4.f

%making a heatmap to compare aerobic and anaerobic conditions
xvalues = {'aerobic','anaerobic'};
yvalues = {'ACALD','ACALDt','ACKr','ACONTa','ACONTb','ACt2r','ADK1','AKGDH','AKGt2r','ALCD2x','ATPM','ATPS4r','Biomass_Ecoli_core_N(w/GAM)-Nmet2','CO2t','CS','CYTBD','D-LACt2',...                        
'ENO','ETOHt2r','EX_ac(e)','EX_acald(e)','EX_akg(e)','EX_co2(e)','EX_etoh(e)','EX_for(e)','EX_fru(e)','EX_fum(e)','EX_glc(e)','EX_gln-L(e)','EX_glu-L(e)','EX_h2o(e)',...
'EX_h(e)','EX_lac-D(e)','EX_mal-L(e)','EX_nh4(e)','EX_o2(e)','EX_pi(e)','EX_pyr(e)','EX_succ(e)','FBA','FBP','FORt2','FORti','FRD7','FRUpts2','FUM','FUMt2_2','G6PDH2r',...
'GAPD','GLCpts','GLNS','GLNabc','GLUDy','GLUN','GLUSy','GLUt2r','GND','H2Ot','ICDHyr','ICL','LDH_D','MALS','MALt2_2','MDH','ME1','ME2','NADH16','NADTRHD','NH4t','O2t','PDH','PFK','PFL','PGI',...
'PGK','PGL','PGM','PIt2r','PPC','PPCK','PPS','PTAr','PYK','PYRt2r','RPE','RPI','SUCCt2_2','SUCCt3','SUCDi','SUCOAS','TALA','THD2','TKT1','TKT2','TPI'};
figure(3)
h=heatmap(xvalues,yvalues,fluxData2);
h.Title = 'Comparing aerobic and anaerobic flux solution vectors';
