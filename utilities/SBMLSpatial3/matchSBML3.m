function [requestedModels,param] = matchSBML3(sbmlFile,modelsdir,param)
%
%This function takes an SBML .xml file and determines which models it
%should match them with. It then returns a list of requested models to be
%synthesized. This is the "cheating" way. 
%this is the 

%load the SBML file
try
    modelObj = sbmlimport(sbmlFile);
catch
    warning('sbmlimport/SimBiology toolbox not available. Parsing SBML using SBMLcompartments3.m')
    modelObj = SBMLcompartments3(sbmlFile);
end

%get the available models in the model directory 
availModels = ml_ls(modelsdir);

nmodels = 1;

reqmodels = {};
j = 1;
Pcount = 0;
Ncount = 0;
Ccount = 0;
supported_cases = ones(1,13);

for i = 1:length(modelObj.Compartments)
    switch upper(modelObj.Compartments(i).name)
        case{'EC','EXTRACELLULAR'}
            %"Extra Cellular" - don't need a model for this one.
            continue
        case{'PM','PLASMAMEM','CELLMEM','CM','CELLMEMBRANE'}
            %Check that we haven't used this before
            supported_cases(2) = supported_cases(2)-1;
            %If we have, skip it! 
            if supported_cases(2)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"CELL MEMBRANE" - need a cell model 
            reqmodels{j} = 'cell';
            
            %assume that protein models contain this
        case{'CP','CYTO','CYTOPLASM'}
            %Check that we haven't used this before
            supported_cases(3) = supported_cases(3)-1;
            %If we have, skip it! 
            if supported_cases(2)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"CYTOPLASM" - need a cell model
            reqmodels{j} = 'cell';
            Ccount = Ccount+1;
            param.SBML_CName{Ccount} = modelObj.Compartments(i).name;
%             ismem(j) = 0;
        case{'NM','NUCMEM','NUCLEARMEM','NUCLEARMEMBRANE'}
            %Check that we haven't used this before
            supported_cases(4) = supported_cases(4)-1;
            %If we have, skip it! 
            if supported_cases(4)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"NUCLEAR MEMBRANE" - need a nuclear model.
            reqmodels{j} = 'nuc';
            ismem(j) = 1;
        case{'NU','NUCLEUS','NUC'}
            %Check that we haven't used this before
            supported_cases(5) = supported_cases(5)-1;
            %If we have, skip it! 
            if supported_cases(5)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"NUCLEUS" - need a nuclear model.
            reqmodels{j} = 'nuc';
            Ncount = Ncount+1;
            param.SBML_NName{Ncount} = modelObj.Compartments(i).name;
%             ismem(j) = 0;
        case{'EM','ENDOMEM','ENDOMEMBRANE','ENDOSOMEMEMBRANE'}
            %Check that we haven't used this before
            supported_cases(6) = supported_cases(6)-1;
            %If we have, skip it! 
            if supported_cases(6)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"ENDOSOMAL MEMBRANE" - need an endosome model
            reqmodels{j} = 'endosome';
            ismem(j) = 1;
        case{'EN','ENDO','ENDOSOME'}
            %Check that we haven't used this before
            supported_cases(7) = supported_cases(7)-1;
            %If we have, skip it! 
            if supported_cases(7)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"ENDOSOME" - need an endosome model
            reqmodels{j} = 'endosome';
            Pcount = Pcount+1;
            param.SBML_PName{Pcount} = modelObj.Compartments(i).name;
%             ismem(j) = 0;
        case{'LM','LYSOMEM','LYSOMEMBRANE','LYSOSOMEMEMBRANE'}
            %Check that we haven't used this before
            supported_cases(8) = supported_cases(8)-1;
            %If we have, skip it! 
            if supported_cases(8)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"LYSOSOMAL MEMBRANE" - need an lysosomal model
            reqmodels{j} = 'lysosome';
            ismem(j) = 1;
        case{'LY','LYSO','LYSOSOME'}
            %Check that we haven't used this before
            supported_cases(9) = supported_cases(9)-1;
            %If we have, skip it! 
            if supported_cases(9)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"LYSOSOME" - need an lysosome model
            reqmodels{j} = 'lysosome';
            Pcount = Pcount+1;
            param.SBML_PName{Pcount} = modelObj.Compartments(i).name;
%             ismem(j) = 0;
        case{'MM','MITOMEM','MITOMEMBRANE','MITOCHONDRIONMEMBRANE'}
            %Check that we haven't used this before
            supported_cases(10) = supported_cases(10)-1;
            %If we have, skip it! 
            if supported_cases(10)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"MITOCHONDRION MEMBRANE" - need an mito model
            reqmodels{j} = 'mitochondrion';
            ismem(j) = 1;
        case{'MITO','MITOCHONDRION'}
            %Check that we haven't used this before
            supported_cases(11) = supported_cases(11)-1;
            %If we have, skip it! 
            if supported_cases(11)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"MITOCHONDRION" - need an mito model
            reqmodels{j} = 'mitochondrion';
            Pcount = Pcount+1;
            param.SBML_PName{Pcount} = modelObj.Compartments(i).name;
%             ismem(j) = 0;
        case{'NUM','NUCLEOMEM','NUCLEOLARMEM','NUCLEOLARMEMBRANE','NUCLEOLIMEM','NUCLEOLIMEMBRANE'}
            %Check that we haven't used this before
            supported_cases(12) = supported_cases(12)-1;
            %If we have, skip it! 
            if supported_cases(12)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"NUCLEOLAR MEMBRANE" - need an nucleolar model
            reqmodels{j} = 'nucleoli';
            ismem(j) = 1;
        case{'NUCLEOLI'}
            %Check that we haven't used this before
            supported_cases(13) = supported_cases(13)-1;
            %If we have, skip it! 
            if supported_cases(13)<0
                warning(['Attempted to duplicate a compartment!',...
                    ' Skipping compartment. One or more compartment will not be mapped.'])
                continue
            end
            %"NUCLEOLI" - need an nucleolar model
            reqmodels{j} = 'nucleoli';
            Pcount = Pcount+1;
            param.SBML_PName{Pcount} = modelObj.Compartments(i).name;
%             ismem(j) = 0;
        otherwise
            warning('Unrecognized/Unsupported model! the resulting SBML spatial file may not be complete');
    end
    
    
    j = j+1;
end

%eliminate duplicates
reqmodels = unique(reqmodels);


%loop through and grab the protein models we need
k = 1;
requestedModels = {};
for i = 1:length(reqmodels)
    %If currmodel is cell or nuc we'll get them at the end so skip
    if strcmpi(reqmodels{i},'cell')||strcmpi(reqmodels{i},'nuc')
        continue
    end
    %otherwise lets find it
    for j = 1:length(availModels)
        %note we assume models are saved as "model"
        load([modelsdir,filesep,availModels{j}]);
        if strcmpi(model.proteinModel.class,reqmodels{i})
            requestedModels{k} = [modelsdir,filesep,availModels{j}];
            ['Success: Model found for ',reqmodels{i},'. Using ',requestedModels{k},' and continuing to next compartment.']
            k = k+1;
            break
        elseif j==length(availModels)
            warning(['Failure: Appropriate model not found for ',reqmodels{i},'. ignoring compartment'])
        end
    end
end

%figure out if cell and nuc models are required 
% cellreq = sum(cellfun(@sum,strfind(reqmodels,'cell')));
% nucreq = sum(cellfun(@sum,strfind(reqmodels,'nuc')));


if isempty(requestedModels)
    requestedModels = availModels{1};
    param.synthesis = 'framework';
else
    param.synthesis = 'all';
end