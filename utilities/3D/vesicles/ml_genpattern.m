function [protimage,nucimage,cellimage,rgbimg] = ml_genpattern(model,param)
%ML_GENPATTERN Generate images for a protein pattern
%   PROTIMG = ML_GENPATTERN(MODEL) retruns the protein image generated from
%   the [generative model] MODEL.
%   
%   PROTIMG = ML_GENPATTERN(MODEL,PARAM) allows users to customize parameters
%   for synthesis.
%   
%   [PROTIMG,DNAIMG,CELLIMG] = ML_GENPATTERN(...) returns DNA image and cell
%   image as well.
%   
%   See also

%   24-Jul-2006 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('imageSize',[1024 1024],'gentex',0,'loc','all','isshow','yes'));

switch model.name
    case 'vesicle'
        [nucEdge,cellEdge] = ml_gencellcomp(model,param);
    
        protimage = ...
            ml_genprotimg(model.proteinModel,nucEdge,cellEdge,param);
        
        nucimage = nucEdge;
        cellimage = cellEdge;
        
        %map nuclear texture
        if exist('nuctex','var')
        end

        se=strel('disk',4,4);
        cellimage2 = imdilate(cellimage,se);
        if ~exist('nuctex','var')
            nucimage2 = imdilate(nucimage,se);
        else
            nucimage2 = nucimage;
        end
        
        rgbimg=ml_synrgbimg(nucimage2,protimage, ...
            cellimage2,'isc');
        
        if strcmp(param.isshow,'yes')
            imshow(rgbimg,[])
        end
end
