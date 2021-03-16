function [requestedModels] = parseSBML3(sbmlFile,modelsdir)
%
%This function takes an SBML .xml file and determines which models it
%should match them with. It then returns a list of requested models to be
%synthesized. 
%this is the 

%load the SBML file
modelObj = sbmlimport(sbmlFile);

%get the available models in the model directory 
availModels = ml_ls(modelsdir);

nmodels = 1;

%loop through the compartments and decide which models are necessary 
for i = 1:length(modelObj.Compartments)
    currCompartment = modelObj.Compartments(i);
    %loop through the available models and determine if there is a suitable
    %model for the current compartment
    for j = 1:length(availModels)
        
        %load currentmodel
        %NOTE: this assumes that the models are named as "model"
        load([modelsdir,availModels{j}]);
        currsubmodels = fieldnames(model);
        for k = 1:length(currsubmodels)
            %figure out whether the current subfield 
            if ~isempty(strfind(lower(currsubmodels{k}),'model'))||~isempty(strfind(lower(currsubmodels{k}),'shape'))
                %check if the submodels have anything in common with the
                if isfield(model.(currsubmodel{k}),'class')
                    compname = model.(currsubmodel{k}).class;
%                 elseif isfield(model.(currsubmodel{k}),'name')
%                     compname = model.(currsubmodel{k}).name;
                else
                    compname = model.(currsubmodel{k});
                end
                if ~isempty(strfind(compname),currCompartment);
                    requestedModels{nmodels} = modelObj.(currsubmodel);
                    nmodels = nmodels+1;
                end
            end
        end
        
    end
    
end