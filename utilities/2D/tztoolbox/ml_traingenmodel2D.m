
function model = ml_traingenmodel2D( protimages, dnaimages, cellimages, ...
    masks, param )
%ML_TRAINGENMODEL2D Train a 2D generative model.
%   MODEL = ML_TRAINGENMODEL(PROTIMAGES,DNAIMAGES,CELLIMAGES,MASKS)
%   returns a generative model that is trained on the images specified by
%   the [string array]s PROTIMAGES (protein channel), DNAIMAGES (DNA
%   channel), CELLIMAGES (cell channel), MASKIMAGES (region masks).
%   MASKIMAGES is empty if there is no mask.
%
%   MODEL = ML_TRAINGENMODEL(PROTIMAGES,DNAIMAGES,CELLIMAGES,MASKS,PARAM)
%   specifies how to train the model according to the structure PARAM:
%

%   17-Jul-2006 Initial write T. Zhao
%   Copyright (c) 2006-2018 Murphy Lab
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

% March 23, 2012 Changed isshow from true to false
%
% March 23, 2012 Renamed the method to ml_traingenmodel2D to avoid
%                 overloading with 3D methods
%
% May 7, 2013 I. Cao-Berg Modifiy algorithm so that if the number of
% images if less than 10 then the number of components is set to the
% number of images
%
% May 7, 2013 I. Cao-Berg Updated method so that it will save the
% temporary results of the cell codes so that it does not need to
% recalculate them
%
% July 7, 2013 G. Johnson fixed bug NaNs in object intensity and 
% covariance matricies were resulting in NaN model parameter vlaues
%
% January 13, 2018 I. Cao-Berg Cleaned up method

if nargin < 3
    error('3 or 4 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

%icaoberg 5/20/2013
param = get_cellorganizer_default_parameters( 'training', param );

%D. Sullivan 6/5/13 imageSize should not be hard coded, need to fix this.
%also, why do we have both param.display and param.disp?
param = ml_initparam(param, ...
    struct('modelname','vesicle'...
    ,'imageSize',[1024 1024] ...
    ,'disp',0 ...
    ,'threshmeth_dna', 'nih' ...
    ,'threshmeth_prot', 'nih'));

model.name = param.modelname;

if ~isfield(param,'ml_objgaussmix')
    param.ml_objgaussmix = struct([]);
end

%icaoberg - changed isshow to 0
param.ml_objgaussmix = ml_initparam(param.ml_objgaussmix, ...
    struct('filter',fspecial('gaussian',10,10),'mindist',10, ...
    'isshow',0, 'gmm',struct('covartype','spherical',...
    'options',[0 1 0 0 0 0 0 0 0 0 0 0 0 500])));

rmidx = [];

temporary_files_folder = [pwd filesep 'temp'];

%icaoberg 5/20/2013
%if these three flags are set to true that means the user would like to
%save the figure shown during debugging
if param.debug && param.display && param.save_helper_figures
    helper_figures_folder = [temporary_files_folder filesep 'figures'];
    subfolders = { 'nucleus', 'cell', 'protein' };
    for s=1:1:length(subfolders)
        subfolder = [helper_figures_folder filesep subfolders{s} ];
        if ~exist( subfolder )
            mkdir( subfolder )
        end
    end
end
clear subfolders
clear s

procfiles = {};

%D. Sullivan 1/30/14 - adding skip if we are just training the cell
%model.(neuron)
if ~strcmpi(param.train.flag,'cell')
    
for i=1:length(dnaimages)
    procfiles{i} = [temporary_files_folder filesep 'preprocessed' ...
        filesep 'cell' num2str(i) '.mat'];    
    if ~exist(procfiles{i})
        disp(['Processing image: ' num2str(i)]);
        
        if ~isempty(masks)
            if ~isempty(masks{i})
                maskimg = ml_readimage(masks{i});
            else
                maskimg = [];
            end
        else
            maskimg = [];
        end
        
        disp('Preprocessing nuclear image');
        %read image
        tic
        img = ml_readimage( dnaimages{i} );
        if isempty( img )
            disp([ 'Unable to load image. Skipping ' dnaimages{i} '.'] );
            continue
        end
            
        %preprocessing
        procimg = ml_preprocess(double(img),maskimg,'ml','yesbgsub', param.threshmeth_dna);
        nucbodyimg = imfill(procimg>0,'hole');
        
        %edge detection (use tz_imedge)
        %   nucedge = bwperim(nucbodyimg);
        nucbody = ml_findmainobj_bw(nucbodyimg);
        toc

        %icaoberg 5/20/2013
        %this display the raw image and the preprocessed image
        if param.debug && param.display
            %this is done to open a maximized figure
            screen_size = get(0, 'ScreenSize');
            f1 = figure;
            set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            
            %make and label plots
            subplot(1,2,1);
            imshow( img );
            title('Raw image');
            subplot(1,2,2);
            imshow( nucbodyimg );
            title('Found nuclear body');
            mtit( dnaimages{i} )
            
            %save helper figures
            if isfield( param, 'save_helper_figures' ) && ... 
                    param.save_helper_figures
                subfolder = [ helper_figures_folder filesep 'nucleus' ];
                [path, filename, extension] = fileparts( dnaimages{i} );
                print( gcf, '-dtiff', [ subfolder filesep filename extension ]);
            end
            clear subfolder
            
            %icaoberg 23/1/2015
            if isfield( param, 'interactive' ) && ...
                    param.interactive
                
                waitforbuttonpress;
            else
                if ~isfield( param, 'time_to_pause' )
                    param.time_to_pause = 2;
                end
                
                pause( param.time_to_pause );
            end
            
            close
        end
        
        disp('Preprocessing cell image');
        if ~isempty( cellimages )
            tic
            %read cell image
            img = ml_readimage(cellimages{i});
            if isempty( img )
                disp([ 'Unable to load image. Skipping ' cellimages{i} '.'] );
                continue
            end
            
            if ~isempty(maskimg)
                img(maskimg==0) = 0;
            end
            
            if length(unique(img(:))) > 2
                celledge = ml_imedge(img,[],'ce');
                cellbodyimg = imfill(celledge,'hole');
            else
                cellbodyimg = img;
            end
            
            cellbody = ml_findmainobj_bw(cellbodyimg);
            toc
        else
            cellbodyimg = [];
            cellbody = [];
        end
        
        cellbody = ml_findmainobj_bw(cellbodyimg);
        toc
        
        %icaoberg 5/20/2013
        %this display the raw image and the preprocessed image
        if param.debug && param.display
            %this is done to open a maximized figure
            screen_size = get(0, 'ScreenSize');
            f1 = figure;
            set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            
            %make and label plots
            subplot(1,2,1);
            imshow( img );
            title('Raw image');
            subplot(1,2,2);
            imshow( cellbodyimg );
            title('Found cell body');
            mtit( dnaimages{i} )
            
            %save helper figures
            if param.save_helper_figures
                subfolder = [ helper_figures_folder filesep 'cell' ];
                [path, filename, extension] = fileparts( cellimages{i} );
                print( gcf, '-dtiff', [ subfolder filesep filename extension ]);
            end
            clear subfolder
            
            if param.interactive
                waitforbuttonpress;
            else
                pause( param.time_to_pause );
            end
            
            close
        end
        
        disp('Preprocessing protein image');
        tic
        if ~isempty( protimages )
            img = ml_readimage(protimages{i});
            if isempty( img )
                disp([ 'Unable to load image. Skipping ' protimages{i} '.'] );
                continue
            end
            
            %icaoberg 5/20/2013
            if param.debug && param.display
                imshow( img );
                title( dnaimages{i} )
                
                if param.interactive
                    waitforbuttonpress;
                else
                    pause( param.time_to_pause );
                end
            end
            

            procimage = ...
                ml_preprocess(double(img),maskimg,'ml','yesbgsub',param.threshmeth_prot);
        else
            procimage = [];
        end
        toc
        
        %icaoberg 5/20/2013
        %this display the raw image and the preprocessed image
        if param.debug && param.display
            %this is done to open a maximized figure
            screen_size = get(0, 'ScreenSize');
            f1 = figure;
            set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            
            %make and label plots
            subplot(1,2,1);
            imshow( img );
            title('Raw image');
            subplot(1,2,2);
            imshow( procimage );
            title('Found protein pattern');
            mtit( dnaimages{i} )
            
            %save helper figures
            if param.save_helper_figures
                subfolder = [ helper_figures_folder filesep 'protein' ];
                [path, filename, extension] = fileparts( protimages{i} );
                print( gcf, '-dtiff', [ subfolder filesep filename extension ]);
            end
            clear subfolder
            
            if param.interactive
                waitforbuttonpress;
            else
                pause( param.time_to_pause );
            end
            
            close
        end
        
        %save results to disc
        save( procfiles{i}, 'nucbody', 'cellbody', ...
            'nucbodyimg', 'cellbodyimg', 'procimage', 'maskimg', 'img' );
        close all
    else
        disp('Preprocessed results found. Skipping recalculation.');
        load( procfiles{i} );  
    end
    
    %NUCLEAR
    disp('Extracting medial axis')
    if ~exist( [temporary_files_folder filesep 'medial_axis'] )
        mkdir([temporary_files_folder filesep 'medial_axis']);
    end
    
    filename = [temporary_files_folder filesep 'medial_axis' ...
        filesep 'cell_' num2str(i) '.mat'];
    
    if ~exist( filename )
        tic
        %rotate nuclear image
        theta = ml_majorangle(nucbodyimg)*180/pi;
        nucbodyimg = ml_rotate(nucbodyimg,-theta);
        
        %extract medial axis
        [imgaxis, medaxis, width] = ...
            ml_imaxis(uint8(nucbodyimg));
        
        %fit medial axis
        if i == 85
            disp('asdf')
        end
        
        shape = ml_mxs2mxp(medaxis,width);
        shapes{i}.shape = shape;
        save( filename, 'shape' )
        toc
    else
        disp(['Medial axis found for image ' num2str(i) '. Skipping recalculation.']);
        tic
        shape = load( filename );
        shapes{i} = shape;
        clear shape
        toc
    end
end

if ~isfield(param,'nucshapemodel')
    param.nucshapemodel = struct([]);
end

param.nucshapemodel = ml_initparam(param.nucshapemodel, ...
    struct('modelname','mxp','constknots',{{0.5,0.5}}));

%Train the statistical model of the features
try
    model.nuclearShapeModel = ml_trainshapemodel2D( shapes, param.nucshapemodel );
catch err
    disp( 'Unable to train nuclear shape model. Returning empty model.' )
    model = [];
    if param.debug
        getReport( err, 'extended' );
    end
    return
end
end
%icaoberg 5/7/13
if strcmpi( param.train.flag, 'nuclear' )
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CELL SHAPE MODEL
for i=1:length(cellimages)
    disp(['Calculation of cell object statistics on image ' num2str(i)]);
    load(procfiles{i}); 
    
    if ~exist([temporary_files_folder filesep 'cellcodes'] )
        mkdir([temporary_files_folder filesep 'cellcodes']);
    end
    
    filename = [temporary_files_folder filesep 'cellcodes' ...
        filesep 'cell_' num2str(i) '.mat'];
    
    if ~exist( filename )
        tic
        
        load(procfiles{i})
        
        param.image_index = i;
        
        %icaoberg 5/13/2013
        try
            anglestep = param.model.anglestep;
        catch
            anglestep = 1;
        end
        
        param.nuclear_image_filename = dnaimages{i};
        param.cell_image_filename = cellimages{i};
        cellcode = ml_parsecell2D({},cellbody,nucbody,anglestep,size(img),...
            {'da','nucarea','nuccenter','nucmangle','nuchitpts',...
            'nuccontour','nucellhitpts','nucdist','nucelldist','nucecc',...
            'cellarea','cellcenter','cellmangle','cellcontour',...
            'cellhitpts','celldist','cellecc'}, param );
        cellcodes{i}.cellcode = cellcode;
        save( filename, 'cellcode' );
        param = rmfield( param, 'image_index' );
        param = rmfield( param, 'cell_image_filename' );
        param = rmfield( param, 'nuclear_image_filename' );
        toc
    else
        disp(['Cell codes found for image ' num2str(i) '. Skipping recalculation.']);
        cellcodes{i} = load( filename );
    end
end

if ~isfield(param,'cellshapemodel')
    param.cellshapemodel = struct([]);
end

%icaoberg 5/7/2013
if length( dnaimages ) >= 10
    number_of_components = 10;
else
    number_of_components = length( dnaimages );
end

param.cellshapemodel = ml_initparam( param.cellshapemodel, ...
    struct('modelname','rds', ...
    'screencell','no',...
    'ml_rdistpca',struct('ml_estpdf', ...
    struct('name','mvn','transform',struct('funname','_pca', ...
    'param',struct('ncomp', number_of_components))), ...
    'startangle','nuc')));

%hidden parameter
try
    param.cellshapemodel.anglestep = param.model.anglestep;
catch
    param.cellshapemodel.anglestep = 1;
end

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
%icaoberg april 11, 2012
model.cellShapeModel = ml_trainshapemodel2D(cellcodes, ...
    param.cellshapemodel );

%icaoberg 5/7/2013
if strcmpi( param.train.flag, 'framework' )
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(protimages)
    disp(['Cell ' num2str(i) '/' num2str(length(protimages))])
    disp('Learning Gaussian mixtures');
    tic
    load(procfiles{i})
    
    objs = ml_findobjs2D(procimage);
    toc
    
    disp(['Objects found: ' num2str(length(objs))]);
    
    if ~exist( [temporary_files_folder filesep 'protein'] )
        mkdir([temporary_files_folder filesep 'protein']);
    end
    
    %PROTEIN MODEL
    for j=1:length(objs)
        disp(['Object:' num2str(j)]);
        obj = objs{j};
        objintens{i}(j) = sum(obj(:,3));
        if size(obj,1)==1
            mixes{i}{j} = [];
        else
            if ~exist( [temporary_files_folder filesep 'protein' ...
                    filesep 'mixture_' num2str(i) '_' num2str(j) '.mat'] )
                tic
                mixture = ml_objgaussmix2D( obj,param.ml_objgaussmix );
                mixes{i}{j} = mixture;
                save( [temporary_files_folder filesep 'protein' ...
                    filesep 'mixture_' num2str(i) '_' num2str(j) '.mat'], ...
                    'mixture', 'obj' );
                
                if param.ml_objgaussmix.isshow == 1
                    drawnow
                end
                
                toc
            else
                disp(['Preprocessed mixture found for image ' num2str(i) ' object ' num2str(j) '.']);
                tic
                load( [temporary_files_folder filesep 'protein' ...
                    filesep 'mixture_' num2str(i) '_' num2str(j) '.mat'] );
                mixes{i}{j} = mixture;
                toc
            end
        end
    end
end

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
                        objintensity*mixes{i}{j}.priors(k)];
                end
            end
        end
    end

    %remove illegal entries
    rminds = isnan(latents) | isnan(intensities);
    latents(rminds) = [];
    intensities(rminds) = [];
    
    if ~isfield(model,'proteinModel')
        model.proteinModel = struct([]);
    end

    model.proteinModel = ml_initparam( model.proteinModel, ...
        struct('name','obj', ...
        'objectModel',struct('name','gauss', ...
        'covartype',param.ml_objgaussmix.gmm.covartype, ...
        'stat',struct('transform',struct('funname','sqrt'),'name','exp') )));

    %PROTEIN MODEL
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
        struct('funname','ml_mappos2D', ...
        'param',struct('order',2,'scale',1))));

    cellNumber = length(objcofs);
    for cellidx=1:cellNumber
        cellcode = cellcodes{cellidx}.cellcode;
        cof = round(objcofs{cellidx});
        cellBoundary = ml_closecurve(cellcode.cellcontour);
        nucBoundary = ml_closecurve(cellcode.nuccontour);

        % beta = tz_logreg([ones(size(normdists,1),1) normdists]);

        cellEdge = ml_obj2img2D(cellBoundary,param.imageSize);
        nucEdge = ml_obj2img2D(nucBoundary,param.imageSize);
        
        try
            [distcodes,coords,angles] = ...
                ml_celldistcode2D(cellEdge,nucEdge,{cof},1,'all');

            normdists = distcodes(:,1)./sum(abs(distcodes(:,1:2)),2);
            %x = [ones(size(normdists,1),1) normdiscts angles angles.^2];
            x = ml_evalfun([normdists angles], ...
                model.proteinModel.positionModel.transform);

            model.proteinModel.positionModel.beta(:,cellidx) = ...
                ml_logreg(x,distcodes(:,3));
        
        catch
            model.proteinModel.positionModel.beta(:,cellidx) = ...
                nan(size(model.proteinModel.positionModel.beta,1), 1);
        end
    end
end

end