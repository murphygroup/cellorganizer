function make_intermediate_results_from_temp_folder( param )
%MAKE_INTERMEDIATE_RESULTS_FROM_TEMP_FOLDER Helper method that takes files
%from the temp folder to make the intermediate results

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

if ~isfield(param, 'verbose')
    param.verbose = true;
end

if ~isfield(param, 'debug')
    param.debug = true;
end

%temp files directories
temp_files_directory = [ pwd filesep 'temp' ];

if strcmpi(param.dimensionality, '2D')
    preprocessed_files_directory = [ temp_files_directory filesep ...
        'preprocessed' ];
    cellcodes_directory = [ temp_files_directory filesep ...
        'cellcodes' ];
    protein_model_temp_files_directory = [ temp_files_directory filesep ...
        'protein' ];
else
    preprocessed_files_directory = [ temp_files_directory filesep ...
        'preprocessing' ];
    cellcodes_directory = [ temp_files_directory filesep ...
        'cell_shape_eigen' ];
    protein_model_temp_files_directory = [ temp_files_directory filesep ...
        'protein_objects_gaussian' ];
end

intermediate_results_folder = [ pwd filesep 'intermediate_results' ];
if ~exist( intermediate_results_folder )
    mkdir( intermediate_results_folder )
end

%get the number of files
files = dir( [  preprocessed_files_directory filesep '*.mat' ] );

for index=1:1:length( files )
    if param.verbose && param.debug
        disp( ['Building intermediate results for image: ' num2str( index )] );
        tic
    end
    
    file = files(index).name;
    
    try
        temporary_filename = [ preprocessed_files_directory filesep file ];
        cell_param.preprocessed = load( temporary_filename );
    catch
        warning( ['Unable to load: ' temporary_file ] );
    end
    clear temporary_filename
    
    if strcmpi( param.train.flag, 'all' ) || strcmpi( param.train.flag, 'framework' )
        if strcmpi( param.nucleus.type, 'diffeomorphic' )
            diffeomorphic_temp_files = [ temp_files_directory ...
                filesep 'diffeomorphic' ];
            
            diffeomorphic_intermediate_files = [ intermediate_results_folder filesep ...
                    'diffeomorphic' ];
            
            if ~exist( diffeomorphic_intermediate_files )
                mkdir( diffeomorphic_intermediate_files );
            end
            
            try
                copyfile( diffeomorphic_temp_files, diffeomorphic_intermediate_files );
            catch
                disp(['Could not copy ' diffeomorphic_temp_files ' to ' diffeomorphic_intermediate_files]);
            end
            
            try
                rmdir( [ diffeomorphic_intermediate_files filesep 'distances' ], 's' );
            catch
                disp(['Could not remove ' diffeomorphic_intermediate_files filesep 'distances' ])
            end
        else
            
            if strcmpi( param.dimensionality, '2d' )
                temporary_file = [ cellcodes_directory filesep ...
                    'cell_' num2str(index) '.mat' ];
            else
                temporary_file = [ cellcodes_directory filesep ...
                    'cellcodes_' num2str(index) '.mat' ];
            end
            
            try
                cell_param = load( temporary_file );
            catch
                warning( ['Unable to load: ' temporary_file ] );
            end
            clear temporary_file
            
            %nuclear features
            if strcmpi( param.dimensionality, '3d' )
                temporary_file  = [ pwd filesep 'temp' filesep ...
                    'nuclearfeats' filesep 'nuclearfeats' num2str(index) '.mat' ];
                
                try
                    cell_param.nucleus.features = load( temporary_file );
                catch
                    warning( ['Unable to load: ' temporary_file ] );
                end
                clear temporary_file
            end
        end
        
        if strcmpi( param.train.flag, 'all' )
            %icaoberg
            %the file structure of the temporary files from the vesicle model
            %is different between dimensionalities
            if strcmpi( param.dimensionality, '2D' )
                search_string = [ protein_model_temp_files_directory ...
                    filesep 'mixture_' num2str(index) '*.mat' ];
                mixtures_temp_files = dir( search_string );
                
                if ~isempty( mixtures_temp_files )
                    mixtures = {};
                    
                    for index2 = 1:1:length( mixtures_temp_files )
                        mixture_file = [ protein_model_temp_files_directory filesep ...
                            mixtures_temp_files(index2).name ];
                        
                        try
                            mixtures{index2} = load( mixture_file );
                        catch
                            warning( ['Unable to load: ' mixture_file ] );
                        end
                        clear mixture_file
                    end
                end
                
                cell_param.protein.mixtures = mixtures;
            else

                %object gaussians
                temporary_file  = [ protein_model_temp_files_directory ...
                    filesep 'object_gaussians' filesep 'gaussobjs_' ...
                    num2str(index) '.mat' ];

                try
                    cell_param.protein.gaussian_objects = load( temporary_file );
                catch
                    warning( ['Unable to load: ' temporary_file ] );
                end
                clear temporary_file

                %object positions
                temporary_file  = [ protein_model_temp_files_directory ...
                    filesep 'object_positions' filesep num2str(index) ...
                    '_Beta.mat' ];

                try
                    cell_param.protein.object_positions = load( temporary_file );
                catch
                    warning( ['Unable to load: ' temporary_file ] );
                end
                clear temporary_file

                %object sizes
                temporary_file  = [ protein_model_temp_files_directory ...
                    filesep 'object_sizes' filesep 'sizefits_' ...
                    num2str(index) '.mat' ];

                try
                    cell_param.protein.object_sizes = load( temporary_file );
                catch
                    warning( ['Unable to load: ' temporary_file ] );
                end
                clear temporary_file

                %object statistics
                temporary_file  = [ protein_model_temp_files_directory ...
                    filesep 'object_stats' filesep 'gaussobjs_' ...
                    num2str(index) '.mat' ];

                try
                    cell_param.protein.object_stats = load( temporary_file );
                catch
                    warning( ['Unable to load: ' temporary_file ] );
                end
                clear temporary_file

                %original objects
                temporary_file  = [ protein_model_temp_files_directory ...
                    filesep 'original_objects' filesep 'obj' ...
                    num2str(index) '.mat' ];

                try
                    cell_param.protein.original_objects = load( temporary_file );
                catch
                    warning( ['Unable to load: ' temporary_file ] );
                end
                clear temporary_file

                %compartment stats
                temporary_file  = [ pwd filesep 'temp' filesep ...
                    'compartment_stats' filesep ...
                    'comparts' num2str(index) '.mat' ];

                try
                    cell_param.compartments.stats = load( temporary_file );
                catch
                    warning( ['Unable to load: ' temporary_file ] );
                end
                clear temporary_file
            end
        end
    end
    
    output_file = [ intermediate_results_folder filesep ...
        'cell_' num2str(index) '.mat' ];
    save( output_file, 'cell_param' )
    clear cell_param
    
    if param.verbose && param.debug
        toc
    end
end
