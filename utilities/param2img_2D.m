function [nucimg, cellimg, protimg, options] = param2img_2D(param, options)
%PARAM2IMG_2D Generates a synthetic image from a parameterization of the cell

% Author: Gregory Johnson
%
% Copyright (C) 2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

if nargin > 2
    error( 'CellOrganizer: Wrong number of input arguments.' );
end

if ~exist( 'param', 'var' )
    options  = [];
end

if isempty( param )
    warning( 'Input parameters cannot be empty. Exiting method' );
    img = [];
    return
end

if ~isfield( options, 'image_size' )
    options = ml_initparam(options,struct('imageSize',[1024 1024]));
else
    options = ml_initparam(options,struct('imageSize', options.image_size));
    options = rmfield( options, 'image_size' );
end

% generate nucleus image, from ml_gencellcomp2D.m
shape = param.nuc.shape;
nucbody = ml_findobjs2D(imfill(ml_mxp2img(shape),'hole'));
nucbody = nucbody{1};

if options.debug && options.display
    [bound_box,box_size]=ml_boundbox2D(nucbody(:,1:2));
    offset = round((options.imageSize-box_size)/2)-bound_box(1,:);
    nucbody_temp = ml_addrow(nucbody(:,1:2),offset);
    temp_image = ml_obj2img2D(nucbody_temp,options.imageSize);
    imshow( temp_image, [] );
    pause(1)
    clear temp_image nucbody_temp offset bound_box box_size
end

try
    anglestep = options.model.anglestep;
catch
    anglestep = 1;
end

cellcode = param.cell.cellcode;

% generate nucleus, cell, (protein) images
if strcmpi( options.synthesis, 'framework' ) || ...
        strcmpi( options.synthesis, 'all' )
    
    %         if options.verbose
    %             disp( 'Generating cell shape' );
    %         end
    
    curve2 = cellcode.nucellhitpts;
    
    curve1 = cellcode.nuchitpts;
    cellBoundary = ml_showpts_2d(curve2,'ln',1);
    nucBoundary = ml_showpts_2d(curve1,'ln',1);
    
    [cellBoundbox,cellBoxsize] = ml_boundbox2D(cellBoundary);
    offset = -cellBoundbox(1,:)+1;
    
    nucBoundary = ml_addrow(nucBoundary,offset);
    cellBoundary = ml_addrow(cellBoundary,offset);
    
    options.imageSize = cellBoxsize;
    
    cellimg = ml_obj2img2D(cellBoundary,options.imageSize);
    nucimg = ml_obj2img2D(nucBoundary,options.imageSize);
else
    % only nucleus image
    curve1 = cellcode.nuchitpts;
    nucBoundary = ml_showpts_2d(curve1,'ln',1);
    [nucBoundbox,nucBoxsize] = ml_boundbox2D(nucBoundary);
    offset = -nucBoundbox(1,:)+1;
    nucBoundary = ml_addrow(nucBoundary,offset);
    
    cellimg = [];
    
    options.imageSize = nucBoxsize;
    nucimg = ml_obj2img2D(nucBoundary,options.imageSize);
end

%make sure the nuclear edge is not empty
%the nuclear edge can never be empty
if isempty( nucimg ) && isempty( cellimg )
    disp( 'Nuclear image is empty. Returning empty framework.' );
    outres = [];
    return
end

if strcmpi( options.synthesis, 'framework' ) || ...
        strcmpi( options.synthesis, 'all' )
    
    if isempty( cellimg )
        disp( 'Cell image is empty. Returning empty framework.' );
        nucimg = [];
        cellimg = [];
        return
    else
        %icaoberg 02/01/2016
        [ croppedcellimg, box ] = cropImgND( cellimg );
        nucimg = nucimg(box(1):box(2),box(3):box(4),:);
        cellimg = cellimg(box(1):box(2),box(3):box(4),:);
    end
else
    if isempty( nucimg )
        disp( 'Nuclear image is empty. Returning empty framework.' );
        nucimg = [];
        cellimg = [];
    else
        %icaoberg 02/01/2016
        [ croppednucimg, box ] = cropImgND( nucimg );
        nucimg = nucimg(box(1):box(2),box(3):box(4),:);
    end
end

% set the output resolution
% to do

%{
    if ~exist('options','var')
        options = struct([]);
    end

    nucedge = nucimg;
    celledge = cellimg;
    options = ml_initparam(options,struct('imageSize',size(nucedge), ...
        'loc','all','minv',1,'minmeth','redraw','savePDF',false));
%}

protimg = [];
end