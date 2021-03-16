function percellparam( imdna_path,imcell_path,...
    improt_path,immask_path,cellnum, param )
%Parameterizes an image of a cell for use in the CellOrganizer suite.
%
%Inputs
%imdna_path - cell array of image files for dna parameter extraction
%imcell_path - cell array of image files for cell parameter extraction
%                   (optional if param.train.flag = 'nuclear')  
%improt_path - cell array of image files for protein parameter extraction
%                   (optional if param.train.flag = 'nuclear' or 'framework')
%immask_path - cell array of image files for masking a single cell
%                   from a given image (these are always optional and can
%                   be passed in as an array of [])
%param             a structure holding possible parameter options. This
%                  struct should also contain all temporary results
%                  folders.(see set_temp_result_folders.m)

% Author: Devin Sullivan 6/13/13
%
% Copyright (C) 2007-2013 Murphy Lab
% Carnegie Mellon University
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu

% Gregory Johnson 6/17/13
%   Modified for single image
%   Shortened some variable names 
%   Formatted for readability


%First perform the segmentation if it doesn't exist
currentSegfile = [param.preprocessingFolder filesep 'cell'...
    num2str(cellnum) '.mat'];

if ~exist( currentSegfile, 'file' )

    seg_cell_and_nucleus( ...
        imdna_path, ...
        imcell_path, ...
        improt_path, ...
        immask_path, ...
        param.downsampling, ...
        param.display, ...
        param, ...
        cellnum );
else
    if param.verbose
        disp( [ 'Temporary file ' currentSegfile ' found. Skipping ' ...
            'image preprocessing.' ] );
    end
end


%Compute Nuclear image feats for each file.
spfeattmp = tp_nucimgfeat(currentSegfile,...
    'cylsurf', param,cellnum );

%Now get the Cell level parameters
if ~strcmpi(param.train.flag,'nuc')
    cellfit_percell( ...
        imdna_path,  ...
        imcell_path, ...
        improt_path, ...
        immask_path, ...
        param, ...
        cellnum );
end

%Now compute the protein models.
if strcmpi(param.train.flag,'all')
    protfit_percell( ...
        imdna_path,  ...
        imcell_path, ...
        improt_path, ...
        immask_path, ...
        param, ...
        cellnum );

    %Lastly use any preprocessed results to create compartment
    %models
    compartment_percell( ...
        improt_path, ...
        param, ...
        cellnum)
end

