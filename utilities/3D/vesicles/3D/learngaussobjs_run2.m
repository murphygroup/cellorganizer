function learngaussobjs_run2( dnaImagesDirectoryPath, ...
    cellImagesDirectoryPath, ...
    proteinImagesDirectoryPath, savedir, param )
% Fit gaussian distributions to the objects

% Tao Peng
%
% Copyright (C) 2012-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% April 26, 2012 I. Cao-Berg Added new flag to save method to allow saving of
%              large matrices
% July 30, 2012 I. Cao-Berg Added a new parameter minObjSize for filtering 
%              objects by size
% August 1, 2012 R.F.Murphy Add debug code; remove unnecessary ct indix
% Jan 30, 2012 I. Cao-Berg Made declaration of temporary folder platform
%              and operating system independent
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
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

if ~exist('param','var')
    param = struct([]);
end

%icaoberg 07/20/2012
%D. Sullivan 6/4/13, shouldn't be fixing the minObjSize as a pixel count.
%This should be dependent on resolution.
%min possible vesicle diameter is ~150A
%this yeilds a minimum volume of ~1.8*10^-6um^3. (Brouillette et al.
%1982)
minObjset = 1.8*10^-6;%um^3
param = ml_initparam(param, ...
    struct('modelname','vesicle','imageSize',[],'display',0, 'minObjSize', minObjset,...
    'debug',false));
% param = ml_initparam(param, ...
%     struct('modelname','vesicle','imageSize',[],'disp',0, 'minObjSize', 10,...
%     'debug',false));
model.name = param.modelname;

%D. Sullivan 6/4/13 - compute voxel size based on model resolution
%adding resolution dependence to minObjSize
voxelvolume = prod(param.model.resolution);%um^3/voxel
%need to convert minObjSize to voxel size
%don't model subvoxel resolutions.
if minObjset<voxelvolume
    minObjSize = 2;%voxels
else
    minObjSize = ceil(param.minObjSize/voxelvolume);%voxels
end

if ~isfield(param,'ml_objgaussmix')
    param.ml_objgaussmix = struct([]);
end

param.ml_objgaussmix = ml_initparam(param.ml_objgaussmix, ...
    struct('filter',fspecial3('Gaussian',[10,10,3]),'mindist',5, ...
    'isshow',0,'gmm',struct('covartype','full',...
    'options',[0 1 0 0 0 0 0 0 0 0 0 0 0 1000 1 20])));

temporary_files = ml_ls( [savedir filesep '*_gaussobjs.mat'] );
if ~isempty( temporary_files )
    disp(['Skipping calculation of ' protype ' Gaussian objects, intermediate results found']);
else
    mixes = {};
    objintens = {};
    offsets = {};
    %icaoberg 1/30/2012
    if isfield(param,'temporaryFilesDirectory')
        temporaryFilesDirectory = param.temporaryFilesDirectory;
    else
        temporaryFilesDirectory = [ pwd filesep 'temp' filesep 'protein_objects_gaussian' filesep 'original_objects' filesep ];
    end
    
    files = ml_dir([temporaryFilesDirectory filesep 'obj*.mat']);
    
    for i = 1:1:length(files)
        %load(['./temp/protein_objects_gaussian/original_objects/' protype '/obj' num2str(i) '.mat'])
        load([temporaryFilesDirectory filesep files{i}]);
        disp(['Image: ' num2str(i)])
        
        %learn Gaussian mixture model
        ct = 0;
        cellcounter = 1;
        for r = 1:size(objects,1)
            for c = 1:size(objects,2)
                blockobjs = objects{r,c};
                if ~isempty(blockobjs)
                    disp(['Block (' int2str(r) ', ' int2str(c) ')'])
                    for j = 1:length(blockobjs)
                        obj = blockobjs{j};
                        
                        %icaoberg 7/30/2012
                        %D. Sullivan 6/4/13 changed to size(obj,1) instead
                        %of length(obj), more robust.
                        if size(obj,1) > minObjSize
                            try
                                disp( [ 'Object:' num2str(j) ' has size ' num2str(length(obj)) ] );
                                ct = ct + 1;
                                disp( ['Count: ' num2str(ct) ] );
                                if ~exist( [savedir filesep files{i} '_' num2str(i) '_' ...
                                        num2str(r) '_' num2str(c) '_' num2str(j) '.mat'] )
                                    objectintensity = sum(obj(:,end));
                                    offset = blockSize*[r-1,c-1,0];
                                    
                                    mix = ml_objgaussmix(obj,...
                                        [],param.ml_objgaussmix);
                                    %warning('Matrix must be positive definite');
                                    
                                    if param.debug && param.display
                                        [tempobjimg,tempoffset] = ml_obj2img(obj,[],{'3d','og'});
                                        figure(1)
                                        for kkk=1:size(tempobjimg,3)
                                            imshow(tempobjimg(:,:,kkk))
                                            pause(0.1)
                                        end
                                        figure(2)
                                        for kkk=1:mix.ncentres
                                            plot(mix.centres(1),mix.centres(2),'r.');
                                            hold on
                                        end
                                        hold off
                                        pause
                                    end
                                    
                                    disp( ['(Mix, Object intensity)=' ...
                                        '(' num2str(mix.ncentres) ...
                                        ',' num2str(objectintensity) ')' ] );
                                    save([savedir filesep files{i} '_' num2str(i) '_' ...
                                        num2str(r) '_' num2str(c) '_' num2str(j) '.mat'], ...
                                        'mix', 'offset', 'objectintensity' );
                                else
                                    disp('Intermediate results found. Loading from disk.');
                                    load([savedir filesep files{i} '_' num2str(i) '_' ...
                                        num2str(r) '_' num2str(c) '_' num2str(j) '.mat'])
                                end
                                
                                if cellcounter == 1
                                    cellmixes = mix;
                                    cellobjintensity = objectintensity;
                                    cellobjoffset = offset;
                                else
                                    cellmixes(cellcounter) = mix;
                                    cellobjintensity(cellcounter) = objectintensity;
                                    cellobjoffset(cellcounter,:) = offset;
                                end
                                cellcounter = cellcounter +1;
                            catch err
                                disp(['Unable to calculate mixture for object:' num2str(j)]);
                            end
                        end
                    end
                end
            end
        end
        %D. Sullivan 6/4/13 - check that we have found some objects 
        if exist('cellmixes','var') && exist('cellobjintensity','var') && exist('cellobjoffset','var')
            mixes{end+1} = cellmixes;
            objintens{end+1}= cellobjintensity;
            offsets{end+1}= cellobjoffset;
        else
            warning('No allowed objects found for file %s Ignoring file.',files{i});
        end
    end    
    save([savedir filesep 'gaussobjs.mat'], ...
        'mixes','objintens','offsets','-v7.3')
end
end
