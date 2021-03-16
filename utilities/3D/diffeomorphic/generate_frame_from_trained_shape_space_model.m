function [frame] = generate_frame_from_trained_shape_space_model(model, position, options)
% Generate a frame from a trained shape space model. Return value is a
% structure with fields interpolated_image, total_wall_time, and total_cpu_time.

% Author: Taraz Buck (tebuck@cmu.edu)
%
% Copyright (C) 201-2012  Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% October 1, 2012 I. Cao-Berg Modified method so that it accepts SLML
% models rather than reading from existing files on disk
% November 13, 2012 I. Cao-Berg Made convergence_registration_error_scale an option
% to the user. The default value is 5E-3.
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

%9/22/13 grj - fixed misspelled variable name

% Not options for the user:
clear default_options

convergence_registration_error_scale = 5e-3;

% Spacing is now defined based on Jieyue's/Aabid's preprocessing
% method (which creates images that are 256 x 256 in the XY plane):
spacing = 1;
% Parameters for windowed LDDMM:
filter_radius = 32 / spacing;
window_diameter = 64;
default_options.filter_radius = filter_radius;
default_options.window_radius = window_diameter;
default_options.kernel_z_radius = filter_radius * spacing / 32;
default_options.maximum_deformation_per_step = [1, 1, 0] .* (2. / spacing) + [0, 0, 1] .* .5;

% Allow option of early stopping during interpolations:
default_options.use_known_distances = false;


if ~exist('options', 'var')
    options = default_options;
else
    option_names = fieldnames(default_options);
    for index = 1:length(option_names)
        current_option = option_names{index};
        if ~isfield(options, current_option)
            options = setfield(options, current_option, getfield(default_options, current_option));
        end
    end
end


% options
% registered_shapes = model.nuclearShapeModel.registered_shapes; %commented out by grj on 9/22/13
reference_index = [];
y2 = [];
convex_hull = [];
tes = [];
distances = [];


%     function [shape] = get_aligned_shape(image_index)
% %         shape = CompressLib.decompress(registered_shapes{image_index}); %commented out by grj on 9/22/13
%         shape = model.nuclearShapeModel.imfunc(image_index)
%     end
% 
%     function [shape] = get_aligned_reference_shape(image_index)
% %         shape = CompressLib.decompress(registered_shapes{reference_index}); %commented out by grj on 9/22/13
%         shape = model.nuclearShapeModel.imfunc(reference_index)
%     end



    function [convergence_error] = get_convergence_error(source, target)
        [M, N, P] = size(source);
        convergence_error = ...
            (std(source(:)) + std(target(:)) * .5)^2 * ...
            (M * N * P) * convergence_registration_error_scale;
    end


registration_options = options;
% Remvoe options specific to this function:
registration_options = rmfield(registration_options, 'use_known_distances');
registration_options.registration_error_function = @get_convergence_error;

if ~isfield(registration_options, 'window_size') & isfield(model.nuclearShapeModel.shape_space_options, 'window_size')
    registration_options.window_size = model.nuclearShapeModel.shape_space_options.window_size;
end


%icaoberg 20/02/2013
%shape_space = {y2, convex_hull, tes};
y2 = model.nuclearShapeModel.positions;
%D. Sullivan 12/14/14 - need to catch if the model does not have the convex
%hull precomputed
if ~isfield(model.nuclearShapeModel,'convex_hull')
    model.nuclearShapeModel.convex_hull =    convhulln( y2 ); %result.convex_hull;
end
convex_hull = model.nuclearShapeModel.convex_hull;
tes = model.nuclearShapeModel.tessellation;

shape_space = { y2, convex_hull, tes};
if options.use_known_distances
    shape_space{end + 1} = distances;
end

frame = render_point_3d_windowed(...
    position, shape_space, ...
    model.nuclearShapeModel.imfunc, 0, registration_options);

% frame = frame.interpolated_image;


end

