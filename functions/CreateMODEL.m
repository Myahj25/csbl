

function model = CreateMODEL(dataTable,prepData)
   
    if width(dataTable) < 2
        error('Input data table must have at least two columns: genes and expression levels.');
    end

    % Create the data structure for the model
    dataStruct.genes = dataTable{:, 1}; % gene names
    dataStruct.tissues = dataTable.Properties.VariableNames(2);
    dataStruct.levels = dataTable{:, 3}; % gene TPM values
    dataStruct.threshold = dataTable{:,4};

    % Create the model using ftINIT
    model = ftINIT(prepData, dataStruct.tissues{1},[], [],dataStruct, ...
        {}, getHumanGEMINITSteps('1+0'), false,true);
end
