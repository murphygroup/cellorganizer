function objectMovementModel(model,param)
%OBJECTMOVEMENTMODEL creates random walks for each object within a cell ignoring conflicts with other objects and cell/nuclear boundaries.
% This is intended to be used when a shape space walk is generated with objects
%
%Inputs:
% model = the model from which we are synthesizing
% param = struct array containing information about where to load shape space walk from and save object motion walks to.
%
%Outputs:
% generates and saves frames of object movements

%Author: Devin Sullivan August 2013
% Copyright (C) 2007-2013  Murphy Lab
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


%try to find frame folder
param = ml_initparam(param,struct('framefolder',['.',filesep,'frames',filesep]));

%assign method for the object movement 
if ~isfield(param,'objmovemethod')
    objparam.method = 'brownian';
else
    objparam.method = param.objmovemethod;
end

%load the original objects
if ~isfield(param,'tempfolder')
    param.tempfolder = './temp';
end
load([param.tempfolder filesep 'OriginalObjects.mat'],...
        'GaussObjects','pos');
%define diffusion and time step from model
% diffusion constant
% Dc = 10e-3;
Dc = model.proteinModel.motion.Dc;

%time step
% dt = 1;
dt = model.proteinModel.motion.dt;

    
%find segmented cells and nuclei
%These folders should not be hard coded in the future, but since it's not
%currently a user accessible option to set them, it's ok. 
donecells = ml_ls([param.framefolder filesep 'Cellwalk' filesep '*.tif']);
donenucs = ml_ls([param.framefolder filesep 'Nucwalk' filesep '*.tif']);
%create a blank image 
tmpimg = ml_readimage(donecells{1});
finalsize_x = floor(param.resolution.cell(1)./param.resolution.objects(1).*size(tmpimg,1));
finalsize_y = floor(param.resolution.cell(2)./param.resolution.objects(2).*size(tmpimg,2));
finalsize_z = floor(param.resolution.cell(3)./param.resolution.objects(3).*size(tmpimg,3));
blankimg = zeros(finalsize_x,finalsize_y,finalsize_z);
% blankimg = ml_readimage(donecells{1}).*0;
% blankimg = AdjustResolutions(blankimg,param.resolution.cell,param.resolution.objects);



%initialize newpositions
newpos = pos;
clear pos

%make output directory
savedir = [param.framefolder filesep 'Protwalk' filesep];
if ~exist( savedir )
    mkdir( savedir )
end

%for each step, move the objects according to the motion model.
for frame = 1:param.walksteps
    %first move all the objects 
    newpos = moveObjects(newpos,Dc,dt,objparam);
    
    %put objects into image
    protimg = ml_imaddobj2(blankimg,GaussObjects,...
        struct('method','replace','pos',newpos,'objectmethod',param.sampling.method));
    
    %Load current cell and nucleus
    currcell = ml_readimage(donecells{frame});
    currcell = AdjustResolutions(currcell,param.resolution.cell,param.resolution.objects);
    currnuc = ml_readimage(donenucs{frame});
    currnuc = AdjustResolutions(currnuc,param.resolution.cell,param.resolution.objects);
    
    %mask image based on the allowed compartment
    %eliminate samples outside the cell.
    switch model.proteinModel.cytonuclearflag
        case {'cyt','cyto'}
            codemask = double(currcell)-double(currnuc);
        case {'nuc','nuclear'}
            codemask = double(currnuc - bwperim(currnuc));
        case 'all'
            codemask = double(currcell - bwperim(currcell));
        otherwise
            error(['Unrecognized location name: ' loc]);
    end

    protimg = protimg.*logical(codemask);
    
    %Save the frame
    img2tif(protimg,[savedir filesep 'frame' num2str(frame) '.tif']);
end