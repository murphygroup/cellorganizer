function param = set_temp_result_folders(param)
%SET_TEMP_RESULT_FOLDERS checks the param struct for user defined temp directories, if not
%specified defaults them 
%
%Author: Devin Sullivan 6/13/13
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


if ~isfield(param,'tempparent')
    param.tempparent = [pwd filesep 'temp'];
end
if ~exist( param.tempparent,'dir' )
    mkdir( param.tempparent );
end


%Preprocessed results (segmented cells and nuclei)
if ~isfield(param,'preprocessingFolder')
    param.preprocessingFolder = [ param.tempparent filesep 'preprocessing' ];
end
if ~exist( param.preprocessingFolder,'dir' )
    mkdir( param.preprocessingFolder );
end

%Per-cell nuclear features
if ~isfield(param,'nuctemppath')
    param.nuctemppath = [param.tempparent filesep 'nuclearfeats'];
end

%Per-cell cell features 
if ~isfield(param,'celltemppath')
    param.celltemppath = [param.tempparent filesep 'cell_shape_eigen'];
end
if ~exist(param.celltemppath,'dir')
    mkdir(param.celltemppath)
end

%Per-cell protein features 
%default parent directory
temporaryProtFolder = [ param.tempparent filesep 'protein_objects_gaussian' ];

%preprocessed objects
if ~isfield(param,'objtemppath')
    param.objtemppath = [ temporaryProtFolder filesep 'original_objects'];
end

%fit objects
if ~isfield(param,'savefitdir')
    param.savefitdir = [ temporaryProtFolder filesep 'object_gaussians'];
end
if ~exist(param.savefitdir,'dir')
    mkdir(param.savefitdir);
end

%object stats
if ~isfield(param,'objstatsdir')
    param.objstatsdir = [temporaryProtFolder filesep 'object_stats' filesep];
end

%object sizes 
if ~isfield(param,'objsizedir')
    param.objsizedir = [temporaryProtFolder filesep 'object_sizes' filesep];
end

%object positions 
if ~isfield(param,'objposdir')
    param.objposdir = [temporaryProtFolder filesep 'object_positions' filesep];
end

if ~isfield(param,'compartmentdir')
    param.compartmentdir = [ param.tempparent filesep 'compartment_stats' filesep];
end
if ~exist(param.compartmentdir,'dir')
    mkdir(param.compartmentdir);
end

