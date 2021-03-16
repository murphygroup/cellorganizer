function [ nucimg, cellimg ] = ml_gencellcomp2D( model, param )
%ML_GENCELLCOMP Generate nuclear edge and cell edge from a model.
% NUCIMG = ML_GENCELLCOMP(MODEL) returns the nuclear edge image generated
% from the generative model MODEL.
%
% NUCIMG = ML_GENCELLCOMP(MODEL,PARAM) specifies how to generate the image.
% PARAM has the following fields:
%     'imageSize' - image size
%
% [NUCIMG,CELLIMG] = ML_GENCELLCOMP(...) also returns cell edge image.
%
% See also

% Ting Zhao
%
% Copyright (c) 2007-2015 Murphy Lab
% Carnegie Mellon University
%
% May 13, 2013 I. Cao-Berg
%  -The method now exploits the synthesis flag, that is, if synthesis is
%   set to nucleus it will only generate and return a nuclear image. If
%   synthesis flag is set to framework or all, then it will synthesize the
%   cell framework, that is the nuclear and cell edge images
%  -The method has been updated so that it will try generating a framework
%  a number of times before failing and returning empty images
% I. Cao-Berg 3/5/2015 Included debug, display and verbose flags
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

%maximum number of times it will try to generate a framework before
%returning empty images

maximum_number_of_tries = 1000;

nucimg = []; cellimg = [];

if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

if ~isfield( param, 'imageSize' )
    param = ml_initparam(param,struct('imageSize',[1024 1024]));
end

param = ml_initparam( param, ...
    struct('debug' , false) );
param = ml_initparam( param, ...
    struct('display' , false) );
param = ml_initparam( param, ...
    struct('verbose' , true) );

%icaoberg 7/1/2013
if isfield( model, 'nuclearShapeModel' )
    shape = ml_genshape2D( model.nuclearShapeModel);
else
    if param.debug
        warning( 'Model did not contain a nuclear shape model. Unable to synthesize framework.' )
    end
    return
end

%Convert the contour to a solid object
nucbody = ml_findobjs2D(imfill(ml_mxp2img(shape),'hole'));
nucbody = nucbody{1};

%icaoberg 3/5/2015
%this helper block will open a display and plot the solid nuclear boundary
if param.debug && param.display
    [bound_box,box_size]=ml_boundbox2D(nucbody(:,1:2));
    offset = round((param.imageSize-box_size)/2)-bound_box(1,:);
    nucbody_temp = ml_addrow(nucbody(:,1:2),offset);
    temp_image = ml_obj2img2D(nucbody_temp,param.imageSize);
    imshow( temp_image, [] );
    pause(1)
    clear temp_image nucbody_temp offset bound_box box_size
end

%The cell code base on the synthesized nucleus
%This is necessary for further steps

%icaoberg 7/1/2013
if ~isfield( model, 'cellShapeModel' )
    if param.debug
        warning( 'Model did not contain a cell shape model. Unable to synthesize framework.' )
    end
    return
end

%icaoberg 5/17/2013
try
    anglestep = model.cellShapeModel.anglestep;
catch
    anglestep = 1;
end

cellcode = struct([]);
cellcode = ml_parsecell2D(cellcode,nucbody,nucbody, anglestep , ...
    param.imageSize,...
    {'da','nuchitpts','nucdist','nuccenter','nucmangle',...
    'nucellhitpts','celldist'}, param);

if strcmpi( param.synthesis, 'framework' ) || ...
        strcmpi( param.synthesis, 'all' )
    
    if param.verbose
        disp( 'Generating cell shape' );
    end
    
    cellcode2 = ml_gencellshape(model.cellShapeModel,cellcode);
    
    curve2 = cellcode2.nucellhitpts;
    
    curve1 = cellcode2.nuchitpts;
    cellBoundary = ml_showpts_2d(curve2,'ln',1);
    nucBoundary = ml_showpts_2d(curve1,'ln',1);
    
    [cellBoundbox,cellBoxsize] = ml_boundbox2D(cellBoundary);
    offset = -cellBoundbox(1,:)+1;
    
    nucBoundary = ml_addrow(nucBoundary,offset);
    cellBoundary = ml_addrow(cellBoundary,offset);
    
    
    param.imageSize = cellBoxsize;
    
    cellimg = ml_obj2img2D(cellBoundary,param.imageSize);
    nucimg = ml_obj2img2D(nucBoundary,param.imageSize);
else
    curve1 = cellcode.nuchitpts;
    nucBoundary = ml_showpts_2d(curve1,'ln',1);
    [nucBoundbox,nucBoxsize] = ml_boundbox2D(nucBoundary);
    offset = -nucBoundbox(1,:)+1;
    nucBoundary = ml_addrow(nucBoundary,offset);
    
    cellimg = [];
    param.imageSize = nucBoxsize;
    nucimg = ml_obj2img2D(nucBoundary,param.imageSize);
end
