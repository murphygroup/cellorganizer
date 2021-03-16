function [ nuc, cell, improt ] = intermed2img(varargin)
%INTERMED2IMG Displays per-cell parameterizations from the intermediate results
%directory

% Gregory Johnson
%
% Copyright (C) 2013-2016 Murphy Lab
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

% Feb. 2, 2016 I. Cao-Berg Updated method to match the refactoring of
% CellOrganizer v2.5
%
% March 4, 2016 I. Cao-Berg Updated method to ignore the synthesis of the
% protein pattern unless param.synthesis is set to 'all'

nuc = [];
cell = [];
improt = [];
img = [];

intermed_dir = varargin{1};

try
    filename = [ intermed_dir filesep 'cell.mat'];
    load( filename )
catch
    warning( ['Unable to load temporary file: ' filename ]);
    return
end
clear filename

try
    filename = [ intermed_dir filesep 'nuc.mat'];
    load( filename )
catch
    warning( ['Unable to load temporary file: ' filename ]);
    return
end
clear filename

%icaoberg 12/2/2013
%just in case we lack some of the parameters
try
    coefs = nuc.spfeat.coefs;
    
    u = (coefs(:,1) + coefs(:,end))/2;
    coefs(:,1) = u;
    coefs(:,end) = [];
    
    model = struct;
    model.dimensionality = '3D';
    model.type = 'framework';
    
    model.nuclearShapeModel = struct;
    model.nuclearShapeModel.type = 'cylindrical_surface';
    model.nuclearShapeModel.height.stat.mu = nuc.spfeat.height;
    model.nuclearShapeModel.height.stat.sigma = 0;
    model.nuclearShapeModel.surface.stat.mu = coefs(:)';
    model.nuclearShapeModel.surface.stat.sigma = 0;
    model.nuclearShapeModel.resolution = [1,1,1];
    
    model.nuclearShapeModel.surface.number = nuc.spfeat.number;
    model.nuclearShapeModel.surface.form = nuc.spfeat.form;
    model.nuclearShapeModel.surface.order = nuc.spfeat.order;
    model.nuclearShapeModel.surface.constknot_phi = [0.2500 0.3750 0.5000 0.6250 0.7500];
    model.nuclearShapeModel.surface.constknot_h = 0.5000;
    
    model.cellShapeModel = struct;
    model.cellShapeModel.type = 'ratio';
    model.cellShapeModel.meanShape.const = cell.rad_ratio;
    model.cellShapeModel.modeShape.const = 0;
    
    param.generatemeanshape = true;
    param.synthesis = 'framework';
    
    [nuc, cell] = model2framework(model, param);
catch
    warning('Unable to synthesize framework.');
    return
end

nuc = flipdim(flipdim(permute(nuc, [2,1,3]), 1),2);
cell = flipdim(flipdim(permute(cell, [2,1,3]), 1),2);

%icaoberg 12/2/2013
%just in case we cannot synthesize protein instances
if strcmpi( param.synthesis, 'all' )
    try
        if isfield(cell_param, 'protein')
            centers = vertcat(cell_param.protein.gaussian_objects.mixes{:}.centres);
            covars = cat(3,cell_param.protein.gaussian_objects.mixes{:}.covars);
            objinten = [cell_param.protein.gaussian_objects.mixes{:}.priors];
            
            imsize = size(cell_param.compartments.stats.indexedimg);
            
            improt = gauss2img(centers, covars, objinten, imsize);
        else
            improt = [];
        end
    catch
        warning( 'Unable synthesize protein pattern image' );
        improt = [];
    end
end

% icaoberg 12/2/2013
% i was unable to find a file containing the compartments field. where is
% this setup?
% img = cell_param.compartments.stats.indexedimg;
% img(img == 3) = 1;
end