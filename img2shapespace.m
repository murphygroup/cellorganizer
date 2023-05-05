function answer = img2shapespace(varargin)
% img2shapespace.m
%
% Train 3D generative SPHARM-RPDM cell shape model and creates a report
% using the trained model.
%
% Input
% -----
% * a directory of raw or synthetic cell shape images (cell shapes file
%       format 'LAM_cell[1-9]_ch1_t1.tif') or OME.TIFFs
% * an options structure to override default options
%       see http://www.cellorganizer.org/docs/2.9.0/chapters/cellorganizer_options.html#img2slml
%       important options include
%       .model.name for name of the model itself
%       .model.filename for name of the file in which to store the model
%       .model.resolution for setting image/model pixel size
%       .downsampling for setting downsampling before model building
%       .shape_evolution for setting whether report includes example
%           shape evolution
%
% Output
% ------
% * a valid SLML model file
% * a report with an embedded shape space

% Soham Chakraborti (sohamc@andrew.cmu.edu)
% R.F. Murphy (murphy@cmu.edu)
%
% Copyright (C) 2020 Murphy Lab
% Computational Biology Department
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
% For additional information visit http://murphylab.cbd.cmu.edu or
% send email to murphy@cmu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isdeployed()
    filename = is_deployed(varargin{1});
    load(filename);
else
    if nargin == 2
        dnaImagesDirectoryPath = varargin{1};
        cellImagesDirectoryPath = varargin{2};
        options = struct([]);
    else
        dnaImagesDirectoryPath = varargin{1};
        cellImagesDirectoryPath = varargin{2};
        options = varargin{3};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% img2slml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultoptions = struct();
defaultoptions.verbose = false;
defaultoptions.debug = false;
defaultoptions.display = false;
defaultoptions.model.id = uuidgen();
%is this needed?  what does it do?
defaultoptions.train = struct( 'flag', 'cell' );
defaultoptions.cell.class = 'cell_membrane';
defaultoptions.cell.type = 'spharm_rpdm';
defaultoptions.labels = 'unique';

options = ml_initparam( options, defaultoptions);

if ~isfield(options,'model.filename')
    options.model.filename = [options.model.name '.mat'];
end    

if ~isfield(options,'model.resolution')
    disp ('No resolution specified for the input images;')
    disp ('Using default of [0.1, 0.1, 0.1]');
    % this is the resolution of the input image
    options.model.resolution = [0.1, 0.1, 0.1];
end

if ~isfield(options,'downsampling')
    disp ('No downsampling specified prior to model creation;')
    disp ('Using default of no downsampling');
    options.downsampling = [1, 1, 1];
end

% postprocess of parameterization: alignment
options.model.spharm_rpdm.postprocess = ~false;
% alignment method: 'major_axis' or 'foe'
options.model.spharm_rpdm.alignment_method = 'major_axis';
% plane of rotation: 'xy', 'xz', 'yz' or 'xyz'
options.model.spharm_rpdm.rotation_plane = 'xy';

% degree of the descriptor
options.model.spharm_rpdm.maxDeg = 31;

% latent dimension for the model
options.model.spharm_rpdm.latent_dim = 15;
options.model.spharm_rpdm.segmindnaImagesDirectoryPathfraction = 0.1;

if exist(dnaImagesDirectoryPath,'var') && ~isempty(dnaImagesDirectoryPath)
    options.model.name = 'SPHARM-RPDM-NucCellShapeModel';
    options.dnaImagesDirectoryPathleus.class = 'nuc_membrane';
    options.dnaImagesDirectoryPathleus.type = 'spharm_rpdm';
    options.model.spharm_rpdm.components = {'nuc', 'cell'};
else
    options.model.spharm_rpdm.components = {'cell'};
    options.model.name = 'SPHARM-RPDM-CellShapeModel';
end

tic; answer = img2slml( '3D', dnaImagesDirectoryPath, ...
    cellImagesDirectoryPath, [], options ); toc,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist( [pwd filesep options.model.filename] )
    answer = slml2info( {[pwd filesep options.model.filename]}, options );
    
    if exist( ['./' options.model.name '-report'])
        rmdir(['./' options.model.name '-report'],'s')
    end
    
    try
        movefile( './report' , ['./' options.model.name '-report']);
        web(['./' options.model.name '-report' filesep 'index.html'])
    catch
        disp('Couldn''t rename report folder and open report in browser');
    end
else
    disp ('Failed to find the model file. Exiting application.' );
    answer = false;
end


if isdeployed
    close all
end

end%img2shapespace
