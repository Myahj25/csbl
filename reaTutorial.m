
%loading model
modelName = 'cardiac_mit_glcuptake_atpmax.mat';
modelsDir = getDistributedModelFolder(modelName);
load([modelsDir filesep modelName])

%setting objective function to ATPS4m
objFun = 'ATPS4m';

%closing exchange reactions 
modelalter = modelCardioMito;
exchanges = {'EX_12dgr_m(e)'
 'EX_arachd(e)'
 'EX_co2(e)'
 'EX_coa(e)'
 'EX_crvnc(e)'
 'EX_cys-L(e)'
 'EX_glc(e)'
 'EX_glu-L(e)'
 'EX_gly(e)'
 'EX_glyc(e)'
 'EX_glyc3p(e)'
 'EX_hdca(e)'
 'EX_lac-L(e)'
 'EX_ocdca(e)'
 'EX_ocdcea(e)'
 'EX_ocdcya(e)'
 'EX_ps_m(e)'
 'EX_urea(e)'};
for i = 1:length(exchanges)
    modelalter = changeRxnBounds(modelalter,exchanges{i},0,'l');
end

%selecting carbon source to be fed into the model at at time
allModels = {};
for i = 1:length(exchanges)
    model = modelalter;
 % change bound of the corresponding exchange reaction using 20 units
    model = changeRxnBounds(model, exchanges{i}, -20, 'l');
 
 % new model name
    str = exchanges{i};
    match = {'-','(e)'};
 for j=1:size(match,2)
        str = strrep(str, match{j},{''});
 end
    newModelName = horzcat(str{1},'_model');
 
 % Combine models in a structure
    allModels.(newModelName) = model;
end
 
[essentialRxn4Models, dataStruct] = essentialRxn4MultipleModels(allModels, objFun);
 
%defining essentiality threshold
essentialityRange = [-100,100]; %negative values only represent absent reactions

%identifying  reactions that are essential in at least one model
numModelsPresent = 1;
rxnsOfInterest_1Model = plotEssentialRxns( essentialRxn4Models, essentialityRange, numModelsPresent);

%identifying  reactions that are essential in at least 11 models
numModelsPresent = 11; 
rxnsOfInterest_11Model = plotEssentialRxns( essentialRxn4Models, essentialityRange, numModelsPresent);

%identifying  reactions that are essential in at least 17 models
numModelsPresent = 17;
rxnsOfInterest_17Model = plotEssentialRxns( essentialRxn4Models, essentialityRange, numModelsPresent);

%identifying  reactions that are essential in all models
numModelsPresent = size(essentialRxn4Models(:,2:end),2);
rxnsOfInterest_allModels = plotEssentialRxns( essentialRxn4Models, essentialityRange, numModelsPresent);

%indentifying reactions that are essential and always present in models
essentialityRange = [0,0]; % only reactions

numModelsPresent = 1;
presentRxnsOfInterest_1Model = plotEssentialRxns( essentialRxn4Models, essentialityRange, numModelsPresent);


%indentifying reactions that are essential and always present in models
essentialityRange = [50,100]; % always positive fluxes
numModelsPresent = 1;
RxnsOfInterest = plotEssentialRxns( essentialRxn4Models, essentialityRange, numModelsPresent);

