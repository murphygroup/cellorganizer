function [nucimg, cellimg, param] = model2diffeomorphicInstance( model, param )
%MODEL2DIFFEORMORPHICINSTANCE Helper method that synthesizes a framework
%from a diffeomorphic model

%icaoberg 10/1/2012
%
% 6/25/13 D. Sullivan - moved this to its own file so it can be called
% externally from model2framework.
% 6/25/13 D. Sullivan - added support for param.randomwalk, a flag that
% indicates that a random walk is desired. 
% 9/6/13 added support for the chunk_start method. 
% 9/22/13 grj - Added missing chunk_finish ^^
% 9/22/13 grj - Moved framefolder variable construction outside of
%              conditional part of code
% 4/21/2014 icaoberg fixed typo in method that referenced tessellation

%
% Copyright (C) 2007-2014  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

options = model.cellShapeModel.shape_space_options;

try
    verbose = param.verbose;
    if isa( verbose, 'logical' )
        verbose = false;
    end
catch
    verbose = false;
end

%icaoberg 9/28/2012
try
    debug = param.debug;
    if isa( debug, 'logical' )
        debug = false;
    end
catch
    debug = false;
end

if ~exist('param', 'var')
    param = [];
end

param = ml_initparam(param, struct( ...
    'framefolder', ['.' filesep 'temp' filesep 'frames_' random_string(10, 2) filesep] ...
    ));

framefolder = param.framefolder;

%if the folder doesn't exist yet, make it
if ~isdir(framefolder)
    mkdir(framefolder);
end


try
    %D. Sullivan 6/25/13 added support for random walk
    if ~isfield(param,'randomwalk')||param.randomwalk==0
        %icaoberg 10/1/2012
        %number_shapes = size(y2, 1);
        if ~isfield(param, 'position') || isempty(param.position)
            positions = model.cellShapeModel.positions;
            positions(any(isnan(positions), 2),:) = [];

            number_shapes = size(positions, 1);

            % Generate a random point inside the convex hull of training shapes in the
            % shape space (probably not uniform sampling):
            random_weights = rand(number_shapes, 1);
            random_weights = random_weights ./ sum(random_weights);

            %icaoberg 10/1/2012
            %param.position = random_weights' * y2;

            param.position = random_weights' * positions;

            % Generate the shape associated with that point:
    %         options = struct();
        end
        
        options.use_compression = false;
    
    %D. Sullivan 6/25/13 added random walk code below
    elseif param.randomwalk==1
        
       
        
        %max steps
        if ~isfield(param,'walksteps')
            warning(['No walkstep number specified (param.walksteps).',...
                'synthesizing only 1 image']);
            maxIter = 1;
        else
            maxIter = param.walksteps;
        end
%         maxIter = 50;
        
        %D. Sullivan 7/16/13
        %check if the work has already been done
        %In the future we should skip to the last synthed shape instead of
        %only skipping when totally done. 
        if isdir([framefolder filesep 'Cellwalk'])
            donecells = ml_ls([framefolder filesep 'Cellwalk' filesep '*.tif']);
            donenucs = ml_ls([framefolder filesep 'Nucwalk' filesep '*.tif']);
            if length(donecells)==param.walksteps
%                 for i = 1:length(donecells)
%                     nucimg{i} = ml_readimage(donenucs{i});
%                     cellimg{i} = ml_readimage(donecells{i});

                    nucimg = ml_readimage(donenucs{1});
                    cellimg = ml_readimage(donecells{1});
                    for z = 1:size(cellimg,3)
                        nucimg(:,:,z) = imfill(nucimg(:,:,z),'holes');
                        cellimg(:,:,z) = imfill(cellimg(:,:,z),'holes');
                    end
%                 end
                return
            end
        end
        
        %save some intermediate results (may want to delete these in the
        %final version)
        savepath = [framefolder filesep 'walk'];
%         savepath = '/Users/devinsullivan/MurphyLab/CellShape/walk';
        
        %set random walk paraemters
        %For right now we will hard code all these options while testing
        %type of walk
        if isfield(model,'dynamic')
            try
                walk_type = model.dynamic.type
            catch
                error('No type associated with a dynamic model')
            end
            
        
        elseif isfield(param,'walk_type')
            walk_type = param.walk_type;
        else 
            walk_type = 'brownian';
%             walk_type = 'willmore';
        end
        
        if strcmpi(walk_type,'density')
            energies = model.cellShapeModel.density';
        elseif strcmpi(walk_type,'willmore')
            energies = model.cellShapeModel.Willmore_energy;
        else
            energies = ones(size(model.cellShapeModel.positions,1),1);
        end
        
        
        %D. Sullivan 7/22/13 - changed Dc and dt to be params of the
        %cellshape motion model
        param = ml_initparam(param,struct('Dc', 10e-3));
        param = ml_initparam(param,struct('dt', 1));
        if isfield(model.cellShapeModel,'motion')
            if isfield(model.cellShapeModel.motion,'Dc')
                Dc = model.cellShapeModel.motion.Dc;
            else
                Dc = param.Dc;
                warning(['No diffusion constant found in motion model,'...
                    'defaulting to given value: ',num2str(Dc)]);
%                 Dc = 10e-3;
            end
            if isfield(model.cellShapeModel.motion,'dt')
                dt = model.cellShapeModel.motion.dt;
            else
                dt = param.dt;
                warning(['No time step found in motion model,'...
                    'defaulting to given valu: ',num2str(dt)]);
%                 dt = 1;
            end
        elseif isfield(model,'dynamic')
            Dc = param.Dc;
            dt = param.dt;
        else
            warning(['Shape space walk requested with no motion model.'...
                'defaulting diffusion = 10e-3, time step = 1']);
            % diffusion constant
            Dc = 10e-3;
            
            %time step
            dt = 1;
        end

        
        %Set up the shape_space
        y2 = model.cellShapeModel.positions;
        convex_hull = model.cellShapeModel.convex_hull;
        
        %icaoberg 21/4/2014
        %fixed typo
        tes = model.cellShapeModel.tessellation;
        
        shape_space = { y2, convex_hull, tes};
        
        %make the function call
         param = ml_initparam(param,struct('walk_type', walk_type));
         param = ml_initparam(param,struct('Dc', Dc));
         param = ml_initparam(param,struct('dt', dt));
         param = ml_initparam(param,struct('mmax', maxIter));
         param = ml_initparam(param,struct('save_filename', savepath));
        [result,success] = generate_walk_path(model,shape_space,energies,param);
        
%         [result, success] = generate_walk_path(model,...
%             shape_space, walk_type, Dc, dt, maxIter, savepath,energies);

        param.position = result.walk_path;
        
    end

  
    
    
    %icaoberg 10/1/2012
    %D. Sullivan 6/25/13 - added support for randomwalks 
    frame = cell(1,size(param.position,1));
    nucimg = cell(1,size(param.position,1));
    cellimg = cell(1,size(param.position,1));

    energy = zeros(1,size(param.position,1));
    if size(param.position,1)>1 && isfield(param,'templateJob')
        mkdir walkscripts
    end
    %D. Sullivan 7/22/13 save walk options so individual jobs can load them
    param = ml_initparam(param,struct('tempdir', [param.framefolder filesep 'temp' filesep]));
    
    if ~exist(param.tempdir, 'dir')
        mkdir(param.tempdir)
    end
    
    if size(param.position,1)>1
        disp(['Random walk detected,',...
            'generating and saving frames of the movie to disk. ',...
            'The method will only return the final frame to save memory']);
        
    end
    for i = 1:size(param.position,1)
        framefile = [param.tempdir filesep 'genframe' num2str(i)];
        [startimage, ~, ~, framefile] = chunk_start(framefile);
        
        if startimage
            try
                currframe = i;
                curr_point = param.position(i,:);
%                 %save an file for the current frame to be read by the script
%                 save(framefile,...
%                     'model','options','curr_point','framefolder','currframe');
                frame = makeFrame(model,options,curr_point,framefolder,currframe);
                chunk_finish(param.tempdir,['genframe' num2str(i)])
            catch err
                disp(['Skipping image ' num2str(i) ' due to error']);
                getReport( err, 'extended')
            end
            
            chunk_finish(framefile)
        else
            disp(['Result file found, Skipping image ' num2str(i) ' ']);
        end
    end
    
    %D. Sullivan 6/25/13 temporary conversion back to traditional CO
    %output (not cellarray) so that no code is broken by this addition
    %     if length(cellimg)==1
    %         cellimg = cellimg{1};
    %         nucimg = nucimg{1};
    %     end
    
    %D. Sullivan 7/22/13
    %initialize counter, this is to hold the job open until the first frame is
    %done
    if ~isempty(dir([param.tempdir filesep '*.tmp']))
        warning(['CellOrganizer has detected .tmp files still exist.',... 
            'Either a process is still working or one or more frame',...
            'failed to generate. Returning empty cell and nuclear image.']);
        cellimg = [];
        nucimg = [];
    else    
        cellimg = ml_readimage([framefolder filesep 'Cellwalk' filesep 'frame1.tif']);
        nucimg = ml_readimage([framefolder filesep 'Nucwalk' filesep 'frame1.tif']);
    end
    
catch err
    getReport( err, 'extended')
    nucimg = [];
    cellimg = [];
end
end
