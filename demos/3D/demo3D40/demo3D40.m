function answer = demo3D40(options)
% demo3D40
%
% Train 3D generative framework model from all LAMP2 images in the Murphy Lab 3D HeLa dataset.
%
% Input
% -----
% * a directory of raw or synthetic nucleus images
% * a directory of raw or synthetic cell shape images
% * a directory of raw or synthetic lysosome images
% * the resolution of the images (all images should have the same
%   resolution)
%
% Output
% ------
% * a valid SLML model file

% Xin Lu
%
% Copyright (C) 2017 Murphy Lab
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
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
answer = false;
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D40' );
disp( 'The estimated running time is 12 minutes. Please wait.' );

pattern = 'LAMP2';
dimensionality = '3D';
options.protein.class = 'lysosome';
options.protein.name = 'lamp2';
options.nucleus.name = 'lamp2';
options.cell.model = 'lamp2';

% generic model options
% ---------------------
options.model.name = 'framework';
options.model.id = num2str(now);

% nuclear shape model options
% ---------------------------
% nucleus.type              Holds the nuclear model type. Default is "cylindrical_surface".
options.nucleus.type = 'cylindrical_surface';
options.nucleus.class = 'nuclear_membrane';
options.nucleus.id = num2str(now);

% cell shape model options
% ------------------------
options.cell.type = 'ratio';
options.cell.id = num2str(now);
options.cell.class = 'cell_membrane';

% protein.cytonuclearflag   (optional) Determines whether the protein pattern will be generated in
%                           the cytosolic space ('cyto'), nuclear space ('nuc') or everywhere ('all').
%                           Default is cyto.
if strcmpi( options.protein.class, 'nuc' )
    options.protein.cytonuclearflag = 'nuc';
else
    options.protein.cytonuclearflag = 'cyto';
end

% other options
% -------------
options.verbose = true;
options.debug = false;
options = ml_initparam( options, struct( ...
    'train', struct( 'flag', 'framework' )));

% documentation
% -------------
options.documentation.description='This model has been trained using demo3D40 from CellOrganizer';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
directory = '../../../images/HeLa/3D/ometiff';

if ~exist( directory )
    warning(['Directory ' directory ' does not exist. Exiting demo.']);
    return
else
    files = dir([directory filesep '3D_HeLa_LAMP2*.ome.tif']);
    files = sort_nat({files.name});
    
    for i = 1:10
        masks{i} = @() flipdim( OME_loadchannel( ...
            [directory filesep files{i}],1),3 );
        dna{i} = @() flipdim( OME_loadchannel( ...
            [directory filesep files{i}],2),3 );
        cell{i} = @() flipdim( OME_loadchannel( ...
            [directory filesep files{i}],3),3 );
    end
    
    if isempty( masks ) || isempty( dna ) || isempty( cell )
        warning('One or several list of image files is/are empty. Exiting demo');
        return
    end
    
    options.model.resolution=[0.049, 0.049, 0.2000];
    options.downsampling = [5,5,1];
    options.model.filename = 'model.xml';
    options.masks = masks;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    answer = img2slml( dimensionality, dna, cell, [], options );
end
