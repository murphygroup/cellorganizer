function [mixes, objintens, offsets, options] = learngaussobjs_run3(objects, blockSize, options)
% Fit gaussian distributions to the objects

% Devin Sullivan - adapted from Tao Peng's learngaussobjs_run2.m
%
% Copyright (C) 2012-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
%
% Changes made prior to the 6/13/13 refactoring:
% April 26, 2012 I. Cao-Berg Added new flag to save method to allow saving of
%              large matrices
% July 30, 2012 I. Cao-Berg Added a new parameter minObjSize for filtering 
%              objects by size
% August 1, 2012 R.F.Murphy Add debug code; remove unnecessary ct indix
% Jan 30, 2012 I. Cao-Berg Made declaration of temporary folder platform
%              and operating system independent
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
% June 4, 2013 D. Sullivan removed fixed minObjSize as a pixel count.
%                          Changed to resolution dependent value from 
%                          (Brouillette et al.1982)
%
%%
% June 13, 2013 D. Sullivan Function no longer requires files, now gets all
%                           its information from the preprocessed results
%                           whos folder structure is contained within param
%                           Performs all computations in a per-cell manner
%                           keeping track of which cell it's on using
%                           'currfile' 
%%
% Changes made after to the 6/13/13 refactoring:
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

if ~exist('options','var')
    options = struct([]);
end

%icaoberg 07/20/2012
%D. Sullivan 6/4/13, shouldn't be fixing the minObjSize as a pixel count.
%This should be dependent on resolution.
%min possible vesicle diameter is ~150A
%this yeilds a minimum volume of ~1.8*10^-6um^3. (Brouillette et al.
%1982)
minObjset = 1.8*10^-6;%um^3
options = ml_initparam(options, ...
    struct('modelname','vesicle','imageSize',[],'display',0, 'minObjSize', minObjset,...
    'debug',false));

%D. Sullivan 6/4/13 - compute voxel size based on model resolution
%adding resolution dependence to minObjSize
voxelvolume = prod(options.model.resolution);%um^3/voxel
%need to convert minObjSize to voxel size
%don't model subvoxel resolutions.
if minObjset<voxelvolume
    minObjSize = 2;%voxels
else
    minObjSize = ceil(options.minObjSize/voxelvolume);%voxels
end

if ~isfield(options,'ml_objgaussmix')
    options.ml_objgaussmix = struct([]);
end

%D. Sullivan 11/2/14 - again, should not be hard coding these options as
%they may be resolution dependent. 
options.ml_objgaussmix = ml_initparam(options.ml_objgaussmix, ...
    struct('filter',fspecial3('Gaussian',[10,10,3]),'mindist',5, ...
    'isshow',0,'gmm',struct('covartype','full',...
    'options',[0 1 0 0 0 1e-5 0 0 0 0 0 0 0 1000 1 20])));

% temporary_files = ml_ls( [savedir filesep '*_gaussobjs.mat'] );
mixes = {};
objintens = {};
offsets = {};


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
                       
                        objectintensity = sum(obj(:,end));
                        offset = blockSize.*[r-1,c-1,0];

                        mix = ml_objgaussmix(obj,...
                            [],options.ml_objgaussmix);
                        %warning('Matrix must be positive definite');

                        if options.debug && options.display
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
%                                         pause
                        end

                        disp( ['(Mix, Object intensity)=' ...
                            '(' num2str(mix.ncentres) ...
                            ',' num2str(objectintensity) ')' ] );


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
if exist('cellmixes','var') && exist('cellobjintensity','var') && exist('cellobjoffset','var')
    mixes{end+1} = cellmixes;
    objintens{end+1}= cellobjintensity;
    offsets{end+1}= cellobjoffset;
else
    warning('No allowed objects found. Ignoring file.');
end

end
