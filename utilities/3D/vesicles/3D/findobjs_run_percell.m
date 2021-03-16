function [objects, centers, blockSize, numBlocksX, numBlocksY, options] = findobjs_run_percell(prot,...
    mask,options)
% Extract objects from 3D HeLa images

% Author: Devin Sullivan - adapted from Tao Peng's findobjs_run.m
%
% Copyright (C) 2011-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
%
%Changes from findobjs_run.m:
% June 15, 2012 G. Johnson Changed ml_findobjs loop to accept images of
%         arbritrary x and y dimension
% July 18, 2012 I. Cao-Berg Updated method to ignore single pixel objects
% August 1, 2012 I. Cao-Berg Fixed a bug in the code that would insert an artifact to the 
%                            end of the image that would later be considered by the 
%                            algorithm as an object
% August 2, 2012 D. Sullivan Added masking of protein image using segcell.
%                            this ensures you find no objects outside of
%                            the segmented cell
% Feb 22, 2013 D. Sullivan   Added protein resolution adjustment so that
%                            gaussians are trained on cubic voxels with the
%                            maximum resolution
% Feb 24, 2013 D. Sullivan   Removed Feb 22 change for speed
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
%%%
%%
% June 7-13 2013 D. Sullivan changed to per-cell type computations. Major
%                            changes include passing a single file
%                            path rather than file path and adding a
%                            'currfile' parameter to track what loop
%                            iteration we are currently on in parallel.
%                            Also passing in a tempdirectory
%                            structure as part of param.
%%
%Changes since refactoring:
% Sept 3 2013 D. Sullivan    slight bug fix on check if verbose for
%                            displaying intermediate message to the user.
%
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



options = ml_initparam(options, struct('smoothxy', 10, 'smoothz', 3));

H = fspecial3('Gaussian',[options.smoothxy,options.smoothxy,options.smoothz]);


prot = double(prot);


procimage = ml_preprocess(double(prot),mask,'ml','yesbgsub','nih');


dimobjimg = imfilter(procimage,H,'same');
centers =  ml_imlocalmax(dimobjimg);

%icaoberg 22/4/2014
%reused blocks
blockSize = size(prot);
%blockSize = 100;
numBlocksX = round(size(procimage,1) / blockSize(1));
numBlocksY = round(size(procimage,2) / blockSize(2));

objects = cell(numBlocksX, numBlocksY);
%D. Sullivan - may be able to parallelize this block in the future 
for r = 1:numBlocksX
    for c = 1:numBlocksY
        disp( ['Finding objects in block {' num2str(r) ',' num2str(c) '}' ] );
        if r == numBlocksX
            %icaoberg 1/8/2012
            rEnd = size(procimage,1)-1;
        else
            rEnd = r*blockSize;
        end
        if c == numBlocksY
            %icaoberg 1/8/2012
            cEnd = size(procimage,2)-1;
        else
            cEnd = c*blockSize;
        end

        %icaoberg 17/7/2012
        blockimage = procimage((r-1)*blockSize+1:rEnd,(c-1)*blockSize+1:cEnd,:);
        objs = ml_findobjs(blockimage);
        objSizes = cellfun('size',objs,1);
        objs = objs(objSizes>1);
        objects{r,c} = objs;
    end
end


