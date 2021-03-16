function proteinModel = param2protein_gaussian_model_2D(protparam, cellcodes, imsizes, options)

latents = cellfun(@(x) x.latents, protparam,'UniformOutput', false);
latents = vertcat(latents{:});

intensities = cellfun(@(x) x.intensities, protparam,'UniformOutput', false);
intensities = vertcat(intensities{:});

%remove illegal entries
rminds = isnan(latents) | isnan(intensities);
latents(rminds) = [];
intensities(rminds) = [];
    
proteinModel = [];

proteinModel = ml_initparam(proteinModel, ...
    struct('name','obj', ...
    'objectModel',struct('name','gauss', ...
    'covartype',protparam{1}.options.gmm.covartype, ...
    'stat',struct('transform',struct('funname','sqrt'),'name','exp') )));

%PROTEIN MODEL
proteinModel.objectModel.stat = ...
    ml_estpdf(latents,proteinModel.objectModel.stat);

%Intensity model
fs{1} = struct('funname','ml_div2','nvar',2);
fs{2} = struct('funname','sqrt');
proteinModel.objectModel.relation = ml_makecomfun(fs);
intensdata = ml_evalfun({intensities,latents}, ...
    proteinModel.objectModel.relation);

% proteinModel.objectModel.intensStatModel = ...
%     tz_trainlk(intensdata,struct('method','norm'));

proteinModel.objectModel.intensStatModel = ...
    ml_estpdf(intensdata,struct('name','norm'));

%Get positions
objcofs = cell(1, length(protparam));

for i=1:length(protparam)
    mix = [protparam{i}.mixes{:}];
    objcofs{i} = vertcat(mix.centres);
end
    


%object number statistics
objnums = [];
for i=1:length(objcofs)
    objnums(i,:) = size(objcofs{i},1);
end
proteinModel.objectModel.numStatModel = ...
    ml_estpdf(objnums,struct('name','gamma'));

%Position model
proteinModel.positionModel.name = 'logistic';
proteinModel.positionModel = ml_initparam( ...
    proteinModel.positionModel,struct('transform', ...
    struct('funname','ml_mappos2D', ...
    'param',struct('order',2,'scale',1))));

cellNumber = length(objcofs);
for cellidx=1:cellNumber
    cellcode = cellcodes{cellidx};
    cof = round(objcofs{cellidx});
    cellBoundary = ml_closecurve(cellcode.cellcontour);
    nucBoundary = ml_closecurve(cellcode.nuccontour);

    % beta = tz_logreg([ones(size(normdists,1),1) normdists]);

    cellEdge = ml_obj2img2D(cellBoundary,imsizes(cellidx,:));
    nucEdge = ml_obj2img2D(nucBoundary,imsizes(cellidx,:));

    try
        [distcodes,coords,angles] = ...
            ml_celldistcode2D(cellEdge,nucEdge,{cof},1,'all');

        normdists = distcodes(:,1)./sum(abs(distcodes(:,1:2)),2);
        %x = [ones(size(normdists,1),1) normdiscts angles angles.^2];
        x = ml_evalfun([normdists angles], ...
            proteinModel.positionModel.transform);

        proteinModel.positionModel.beta(:,cellidx) = ...
            ml_logreg(x,distcodes(:,3));

    catch
        proteinModel.positionModel.beta(:,cellidx) = ...
            nan(size(proteinModel.positionModel.beta,1), 1);
    end
end