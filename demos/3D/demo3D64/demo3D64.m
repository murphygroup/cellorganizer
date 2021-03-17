function answer = demo3D64( options )
% demo3D64
%
% Synthesize one 3D image with nuclear, cell shape and a vesicular channel and
% test SBML Spatial export options by saving files for multiple combinations of
% settings.
%
% Input 
% -----
% * a valid CellOrganizer model
%
% Output
% ------
% * 20 SBML instances (`demo3D64/*.xml`)
%   * Filenames indicate a combination of SBML Spatial-related options:
%     * `_SpatialUseAnalyticMeshes` in a filename means
%       `options.output.SBMLSpatialUseAnalyticMeshes = true`
%     * `_SpatialImage` means `options.output.SBMLSpatialImage = true`
%     * `_SpatialVCellCompatible` means
%       `options.output.SBMLSpatialVCellCompatible = true`
%     * `_SpatialUseCompression` means
%       `options.output.SBMLSpatialUseCompression = true`
%     * `_primitives` means `options.output.primitives = true`
%   * Not all options or combinations of options are intended for the end user.
%     We recommend using:
%     * Only `options.output.SBMLSpatialUseAnalyticMeshes = true` if you need
%       geometry in triangle mesh format. We are working on this being
%       CellBlender compatible.
%     * Only `options.output.SBMLSpatialImage = true` if you need geometry in
%       image format.
%   * Not all combinations of options will produce valid SBML output.
%     * `options.output.SBMLSpatialVCellCompatible = true` attempts to format
%       output to be opened in Virtual Cell. Output will not necessarily be
%       valid SBML.
%   * Not all combinations of options will produce valid SBML Spatial output.
%     * `options.output.SBMLSpatialUseAnalyticMeshes = false` uses Matlab's
%       `isosurface`, which does not guarantee manifold, watertight meshes.
%       This sometimes leads to invalid SBML Spatial when there are duplicate
%       faces, for example.
% * single channel TIF files

% Author: Ivan E. Cao-Berg, Taraz Buck
%
% Copyright (C) 2016-2020 Murphy Lab
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

use_profiling = false;
% use_profiling = true;


% Debug flag defaults
debug_recompute_finished_sets = false;
debug_allow_exceptions = false;
debug_options_sets_prefixes = {};
debug_skip_vcell_compatible = false;


% Debug flags
debug_recompute_finished_sets = true;
debug_allow_exceptions = true;
% debug_skip_vcell_compatible = true;


% Print debug flags
debug_recompute_finished_sets
debug_allow_exceptions
debug_options_sets_prefixes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
start_time = tic;
start_cputime = cputime;

if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( mfilename );
disp( 'The estimated running time is about 9 minutes. Please wait.' );

options_base.seed = 12345;
options_base.targetDirectory = pwd;
options_base.numberOfSynthesizedImages = 20;
options_base.numberOfGaussianObjects = 5;
options_base.prefix = 'img';
options_base.compression = 'lzw';
options_base.microscope = 'none';
options_base.sampling.method = 'disc';
options_base.verbose = true;
options_base.debug = false;
options_base.output.tifimages = true;
options_base.output.indexedimage = true;
% options_base.overwrite_synthetic_instances = false;
options_base.rendAtStd = 1.0;
options_base.objstd = options_base.rendAtStd+0.3;
% options_base.overlapsubsize = 1;
% options_base.output.SBML = true;
options_base.output.SBMLSpatial = true;
options_base.output.SBMLSpatialImage = false;
options_base.output.SBMLSpatialUseAnalyticMeshes = false;
% options_base.output.SBMLDownsampling = 1/4;
options_base.output.SBMLDownsampling = [1/16, 1/16, 1];
% options_base.output.SBMLSpatialVCellCompatible = true;
options_base.output.SBMLSpatialUseCompression = false;
options_base.SBML.translations = {'cell', 'CP'; 'nucleus', 'NU'; 'nuc', 'NU'; 'tfr_mat', 'EN'; 'CP_EC', 'PM'; 'CP_EN', 'EM'; 'CP_NU', 'NM'};
options_base.output.primitives = false;
% options_base.output.OMETIFF = true;


% 'cell' and 'nucleus' options unneeded
options_synthesis_choices = {'framework', 'all'};

shape_format_choices = {'analyticMesh', 'image', 'isosurface'};

options_sets = {};

ks = [false, true];
if debug_skip_vcell_compatible
    ks = [false];
end

options = options_base;
for i = 1:length(options_synthesis_choices)
    options.synthesis = options_synthesis_choices{i};
for j = 1:length(shape_format_choices)
    options.output.SBMLSpatialUseAnalyticMeshes = false;
    options.output.SBMLSpatialImage = false;
    options.output.primitives = false;
    switch shape_format_choices{j}
        case 'analyticMesh'
            options.output.SBMLSpatialUseAnalyticMeshes = true;
        case 'image'
            options.output.SBMLSpatialImage = true;
        case 'isosurface'
            % Do nothing
        case 'primitives'
            options.output.primitives = true;
        otherwise
            % Do nothing
    end
for k = ks
    options.output.SBMLSpatialVCellCompatible = k;
for m = [false, true]
    options.output.SBMLSpatialUseCompression = m;
    
    if options.output.primitives && (options.output.SBMLSpatialImage || any(strcmp(options.synthesis, {'cell', 'nucleus', 'framework'})))
        continue;
    end
    
    options.prefix = options.synthesis;
    options_sets{end+1} = options;

end
end
end
end


% Create a unique prefix for each options set
for i = 1:length(options_sets)
    options = options_sets{i};
    options.prefix = options.synthesis;
    if options.output.SBMLSpatialUseAnalyticMeshes
        options.prefix = [options.prefix, '_SpatialUseAnalyticMeshes'];
    end
    if options.output.SBMLSpatialImage
        options.prefix = [options.prefix, '_SpatialImage'];
    end
    if options.output.SBMLSpatialVCellCompatible
        options.prefix = [options.prefix, '_SpatialVCellCompatible'];
    end
    if options.output.SBMLSpatialUseCompression
        options.prefix = [options.prefix, '_SpatialUseCompression'];
    end
    if options.output.primitives
        options.prefix = [options.prefix, '_primitives'];
    end
    options_sets{i} = options;
end


options_sets_prefixes = cellfun(@(x)x.prefix, options_sets, 'UniformOutput', false);
options_sets_prefixes = options_sets_prefixes';
options_sets_prefixes

if length(debug_options_sets_prefixes) > 0
    options_sets2_prefixes = debug_options_sets_prefixes;
else
    options_sets2_prefixes = options_sets_prefixes;
end


if size(options_sets2_prefixes, 2) > 1
    options_sets2_prefixes = options_sets2_prefixes';
end
options_sets2_prefixes

options_sets2 = {};
for i = 1:length(options_sets)
    if any(strcmp(options_sets{i}.prefix, options_sets2_prefixes))
        options_sets2{end+1} = options_sets{i};
    end
end
options_sets = options_sets2;
options_sets_count = length(options_sets)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK

% note, the framework will be synthesized using the first protein model
% found to match the given patterns in the SBML file.
% Changing the order/priority of this is not supported at this time.
list_of_models = {'../../../models/3D/tfr.mat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iteration(given_options_set, given_xml_source_filename, given_xml_destination_filename, given_xml_source_zip_filename, given_xml_destination_zip_filename, given_png_source_filename, given_png_destination_filename)
    if use_profiling
        profile off; profile('-memory','on'); profile off; profile on;
    end
    
    try
        state = rng( given_options_set.seed );
    catch err
        rand( 'seed', given_options_set.seed ); %#ok<RAND>
    end
    
    answer = slml2img( list_of_models, given_options_set );
    
    copyfile(given_xml_source_filename, given_xml_destination_filename);
    copyfile(given_xml_source_zip_filename, given_xml_destination_zip_filename);
    copyfile(given_png_source_filename, given_png_destination_filename);
    
    if use_profiling
        profile_results = profile('info');
        profile_filename = [given_options_set.prefix, '_profile_results', '.mat'];
        profsave(profile_results, [profile_filename, '_html']);
        save(profile_filename, 'profile_results');
    end
    
    %{
    try
        [sbml_valid, sbml_problems] = validate_SBML_instance(xml_destination_filename, true)
        % [sbml_valid, sbml_problems] = validate_SBML_instance(xml_source_zip_filename, true)
    catch ME
        warning('Caught exception:');
        getReport(ME)
    end
    %}
end

for i = 1:length(options_sets)
    options_set = options_sets{i}
    options_set_output = options_set.output
    
    xml_source_filename = [options_set.prefix, '/cell1/cell.xml'];
    xml_destination_filename = [options_set.prefix, '.xml'];
    xml_source_zip_filename = [xml_source_filename, '.zip'];
    xml_destination_zip_filename = [xml_destination_filename, '.zip'];
    png_source_filename = [options_set.prefix, '/cell1/indexed.png'];
    png_destination_filename = [options_set.prefix, '.png'];
    xml_destination_validation_filename = [xml_destination_filename, '.validationresults.xml'];
    % xml_destination_zip_filename = [xml_destination_zip_filename, '.validationresults.xml'];
    
    % if exist(xml_destination_validation_filename, 'file') && ~debug_recompute_finished_sets
    if exist(xml_destination_filename, 'file') && ~debug_recompute_finished_sets
        % Do not repeat work
        continue;
    end
    
    iteration_handle = @()iteration(options_set, xml_source_filename, xml_destination_filename, xml_source_zip_filename, xml_destination_zip_filename, png_source_filename, png_destination_filename);
    
    if debug_allow_exceptions
        iteration_handle();
    else
        try
            iteration_handle();
        catch ME
            warning('Caught exception:');
            getReport(ME)
            % delete(xml_destination_filename);
            % delete(xml_destination_zip_filename);
            if exist(xml_destination_validation_filename, 'file')
                delete(xml_destination_validation_filename);
            end
        end
    end
end

% Create table of features and SBML validity

elapsed_time = toc(start_time);
elapsed_cputime = cputime - start_cputime;
fprintf('\n%s took %.3f s (%.3f s CPU time)\n\n', mfilename, elapsed_time, elapsed_cputime);

end%demo3D64
