function model = ml_traingenmodel(protimages,dnaimages,cellimages, ...
    masks,param)
%ML_TRAINGENMODEL Train a generative model.
%   MODEL = ML_TRAINGENMODEL(PROTIMAGES,DNAIMAGES,CELLIMAGES,MASKS)
%   returns a generative model that is trained on the images specified by
%   the [string array]s PROTIMAGES (protein channel), DNAIMAGES (DNA
%   channel), CELLIMAGES (cell channel), MASKIMAGES (region masks).
%   MASKIMAGES is empty if there is no mask.
%   
%   MODEL = ML_TRAINGENMODEL(PROTIMAGES,DNAIMAGES,CELLIMAGES,MASKS,PARAM) 
%   specifies how to train the model according to the structure PARAM:
%          
%
%   See also

%   17-Jul-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu

if nargin < 3
    error('3 or 4 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('modelname','vesicle','imageSize',[1024 1024],'disp',0));
model.name = param.modelname;

if ~isfield(param,'ml_objgaussmix')
    param.ml_objgaussmix = struct([]);
end

param.ml_objgaussmix = ml_initparam(param.ml_objgaussmix, ...
    struct('filter',fspecial('gaussian',10,10),'mindist',10, ...
    'isshow',1, 'gmm',struct('covartype','spherical',...
    'options',[0 1 0 0 0 0 0 0 0 0 0 0 0 500])));

%train dna shapes   
if param.disp
    disp('train nuclear shapes ...')
    disp('processing images ...')
end

rmidx = [];

for i=1:length(dnaimages)
    if param.disp
        disp(i)
    end
    if ~isempty(masks)
        if ~isempty(masks{i})
            maskimg = ml_readimage(masks{i});
        else
            maskimg = [];
        end
    else
        maskimg = [];
    end
    
    if param.disp
        disp('dna image')
    end
    
    %read image 
    img = ml_readimage(dnaimages{i});

    %preprocessing
    procimg = ml_preprocess(double(img),maskimg,'ml','yesbgsub');
    nucbodyimg = imfill(procimg>0,'hole');
    
    %edge detection (use tz_imedge)
%     nucedge = bwperim(nucbodyimg);
    nucbody = ml_findmainobj_bw(nucbodyimg);
    
    if param.disp
        disp('cell image')
    end
    
    %read cell image
    img = ml_readimage(cellimages{i});
    
    if ~isempty(maskimg)
        img(maskimg==0) = 0;
    end
    
    %preprocessing
    %tz- 06-Oct-2006
%     medimg=medfilt2(img,[5 5]);
%     cellbodyimg = imfill(medimg>0,'hole');
%     celledge = bwperim(cellbodyimg);  
    %tz--
    
    %tz+ 06-Oct-2006
    celledge = ml_imedge(img,[],'ce');
    cellbodyimg = imfill(celledge,'hole');
    %tz++
    
    cellbody = ml_findmainobj_bw(cellbodyimg);
    
    %parse cell
    cellcodes{i} = ml_parsecell({},cellbody,nucbody,1,size(img),...
        {'da','nucarea','nuccenter','nucmangle','nuchitpts',...
        'nuccontour','nucellhitpts','nucdist','nucelldist','nucecc',...
        'cellarea','cellcenter','cellmangle','cellcontour',...
        'cellhitpts','celldist','cellecc'},0);  
    
    if cellcodes{i}.succ==0
        rmidx = [rmidx i];
    end
    
    %rotate nuclear image
    theta = ml_majorangle(nucbodyimg)*180/pi;
    nucbodyimg = ml_rotate(nucbodyimg,-theta);
    

    %extract medial axis
    [imgaxis,medaxis,width] = ...
        ml_imaxis(uint8(nucbodyimg));

    %fit medial axis
    shapes{i} = ml_mxs2mxp(medaxis,width);
    
    if ~isempty(protimages)
        if param.disp
            disp('protein image');
        end
        
        img = ml_readimage(protimages{i});
        
        procimage = ...
            ml_preprocess(double(img),maskimg,'ml','yesbgsub','nih');
        objs = ml_findobjs(procimage);
        
        %learn Gaussian mixture model
        for j=1:length(objs)
            obj = objs{j};
            objintens{i}(j) = sum(obj(:,end));
            if size(obj,1)==1
                mixes{i}{j} = [];
            else
                mixes{i}{j} = ml_objgaussmix(obj,param.ml_objgaussmix);
                if param.ml_objgaussmix.isshow==1
                    drawnow
                end     
            end     
        end  
    end
end

if ~isfield(param,'nucshapemodel')
    param.nucshapemodel = struct([]);
end

param.nucshapemodel = ml_initparam(param.nucshapemodel, ...
    struct('modelname','mxp','constknots',{{0.5,0.5}}));

%train the statistical model of the features
model.nuclearShapeModel = ml_trainshapemodel(shapes,param.nucshapemodel); 

%cell shape model
if ~isfield(param,'cellshapemodel')
    param.cellshapemodel = struct([]);
end

param.cellshapemodel = ml_initparam(param.cellshapemodel, ...
    struct('modelname','rds', 'screencell','no',...
    'ml_rdistpca',struct('ml_estpdf', ...
    struct('name','mvn','transform',struct('funname','_pca', ...
    'param',struct('ncomp',10))), ...
    'startangle','nuc')));

cellIndices = 1:length(dnaimages);

if ~isempty(rmidx)
    cellIndices(rmidx) = [];
    objintens(rmidx) = [];
    mixes(rmidx) = [];
end

combfeats = ml_combccfeats(cellcodes(cellIndices)); 

if strcmp(param.cellshapemodel.screencell,'yes')==1
    %train cell shape
    %Area ratios
    % combratio = combfeats(goodCellIndices,2);
    combratio = combfeats(:,2);

    %Gaussian mixture of area ratios
    [bestk,bestpp,bestmu,bestcov] = ...
        ml_mixtures4(combratio',1,3,0,1e-5,0);

    if bestk>1
%         pp1=bestpp(1);
%         pp2=bestpp(2);
        sigma1=sqrt(bestcov(:,:,1));
        sigma2=sqrt(bestcov(:,:,2));
        mu1=bestmu(1);
        mu2=bestmu(2);

        %Find threshold for small and big ratios
        x = ml_gscross(mu1,sigma1,mu2,sigma2);

        sizeThreshold = x(1);

        smallCellIndices = cellIndices(combratio<sizeThreshold);
%         largeCellIndices = cellIndices(combratio>=sizeThreshold);
    else
        smallCellIndices = cellIndices;
    end
else
    smallCellIndices = cellIndices;
end

%PCA ratio model
%model.cellShapeModel = tz_rdistpca(cellcodes(smallCellIndices),param);

cellcodes = cellcodes(smallCellIndices);

%Cell shape model
model.cellShapeModel = ml_trainshapemodel(cellcodes, ...
    param.cellshapemodel);

if ~isempty(protimages)
    %Parameters of Gaussian objects
    latents = [];
    intensities = [];
    for i=1:length(mixes)
        for j=1:length(mixes{i})
            objintensity = objintens{i}(j);
            if ~isempty(mixes{i}{j})
                for k=1:mixes{i}{j}.ncentres
                    latents = [latents; mixes{i}{j}.covars(1,1,k)];
                    intensities = [intensities; ...
                        objintens{i}(j)*mixes{i}{j}.priors(k)];
                end
            end
        end
    end

    if ~isfield(model,'proteinModel')
        model.proteinModel = struct([]);
    end

    model.proteinModel = ml_initparam(model.proteinModel, ...
        struct('name','obj', ...
        'objectModel',struct('name','gauss', ...
        'covartype',param.ml_objgaussmix.gmm.covartype, ...
        'stat',struct('transform',struct('funname','sqrt'),'name','exp') ...
        )));

    %Protein model
    model.proteinModel.objectModel.stat = ...
        ml_estpdf(latents,model.proteinModel.objectModel.stat);

    %Intensity model
    fs{1} = struct('funname','ml_div2','nvar',2);
    fs{2} = struct('funname','sqrt');
    model.proteinModel.objectModel.relation = ml_makecomfun(fs);
    intensdata = ml_evalfun({intensities,latents}, ...
        model.proteinModel.objectModel.relation);
    % model.proteinModel.objectModel.intensStatModel = ...
    %     tz_trainlk(intensdata,struct('method','norm'));
    model.proteinModel.objectModel.intensStatModel = ...
        ml_estpdf(intensdata,struct('name','norm'));

    %Get positions
    objcofs = {};

    for i=1:length(mixes)
        objcofs{i} = [];
        for j=1:length(mixes{i})
            if ~isempty(mixes{i}{j})
                for k=1:mixes{i}{j}.ncentres
                    objcofs{i}(end+1,:) = mixes{i}{j}.centres(k,:);
                end
            end
        end
    end

    %object number statistics
    objnums = [];
    for i=1:length(objcofs)
        objnums(i,:) = size(objcofs{i},1);
    end
    model.proteinModel.objectModel.numStatModel = ...
        ml_estpdf(objnums,struct('name','gamma'));

    %Position model
    model.proteinModel.positionModel.name = 'logistic';
    model.proteinModel.positionModel = ml_initparam( ...
        model.proteinModel.positionModel,struct('transform', ...
        struct('funname','ml_mappos', ...
        'param',struct('order',2,'scale',1))));

    cellNumber = length(objcofs);
    for cellidx=1:cellNumber
        cellcode = cellcodes{cellidx};
        cof = round(objcofs{cellidx});
        cellBoundary = ml_closecurve(cellcode.cellcontour);
        nucBoundary = ml_closecurve(cellcode.nuccontour);

        % beta = tz_logreg([ones(size(normdists,1),1) normdists]);

        cellEdge = ml_obj2img(cellBoundary,param.imageSize);
        nucEdge = ml_obj2img(nucBoundary,param.imageSize);
        [distcodes,coords,angles] = ...
            ml_celldistcode(cellEdge,nucEdge,{cof},1,'all');
        normdists = distcodes(:,1)./sum(abs(distcodes(:,1:2)),2);
        %x = [ones(size(normdists,1),1) normdists angles angles.^2];
        x = ml_evalfun([normdists angles], ...
            model.proteinModel.positionModel.transform);
        
        model.proteinModel.positionModel.beta(:,cellidx) = ...
            ml_logreg(x,distcodes(:,3));
    end
end
