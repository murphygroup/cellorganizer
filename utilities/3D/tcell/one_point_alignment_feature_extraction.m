function [feature] = one_point_alignment_feature_extraction(align_file, param)
% this function is for extraction features for inference of synapse
% (transformation)
% 
% Author: Xiongtao Ruan
% Date: Aug. 9, 2016


if nargin < 3
    param = struct();
end

param = process_options_structure( ...
        struct('method', 'around_synapse', ...
               'shape', 'sphere', ...
               'radius', 7, ...
               'thickness', 5, ...
               'dimension', 100 ...
               ), param);

if strcmp(param.method, 'around_synapse')
    if ~isfield(param, 'shape_param')
        param.shape_param = struct();
    end
    switch param.shape
        case 'cylinder'
            param.shape_param = process_options_structure(...
                struct('radius', 7, ...
                        'thickness', 5, ...
                        'dimension', 250, ...
                        'sampling_method', 'distance' ...
                        ), param.shape_param);
        case 'sphere'
            param.shape_param = process_options_structure(...
                struct('radius', 6, ...
                        'dimension', 250, ...
                        'sampling_method', 'distance' ...
                        ), param.shape_param);
    end
end

seg_image = align_file.cropped_segmentation_image;
raw_image = align_file.cropped_segmentation_image;
current_landmarks = align_file.current_landmarks;
template_synapse = [36 36 18];
template_centroid = [36 44.3282 18];
template_landmarks = [template_centroid; template_synapse];

switch param.method
    case 'around_synapse'
        feature = one_point_alignment_around_synapse_feature(seg_image, raw_image, template_landmarks, param);
    case 'SURF'
        flag = 1;
        I_2d = reshape_2d(raw_image);
        points = detectSURFFeatures(I_2d);
        [f, vpts] = extractFeatures(I_2d, points);
        
end
        
    
end

function [feature] = one_point_alignment_around_synapse_feature(seg_image, raw_image, landmarks, param)

synapse = landmarks(2, :);
centroid = landmarks(1, :);
shape_param = param.shape_param;

switch param.shape
    case 'cylinder'
        axis = centroid - synapse;
        % convert to right-hand side, xyz -> yxz
        axis_r = axis;
        axis_r(1) = axis(2);
        axis_r(2) = axis(1);
        synapse_r= synapse;
        synapse_r(1) = synapse(2);
        synapse_r(2) = synapse(1);
        T = param.thickness;
        R = 8;
        [cylinder_coordinates] = cylinder_coordinates_extract(synapse_r, axis_r, T, R);
        Xq = cylinder_coordinates(:, 2);
        Yq = cylinder_coordinates(:, 1);
        Zq = cylinder_coordinates(:, 3);
        
        % visualize choosing region
        if false
            close all;
            round_coordinates = round(cylinder_coordinates);
            binary_img_coordinates = zeros(size(seg_image));
            binary_img_coordinates(round_coordinates(:, 1), round_coordinates(:, 2), round_coordinates(:, 3)) = 1;
            
            binary_img_coordinates_2d = reshape_2d(binary_img_coordinates);
            % se = strel('diamond', 2);
            % binary_img_coordinates_2d = imdilate(binary_img_coordinates_2d, se);
            image_to_show = [];
            image_to_show = [image_to_show; reshape_2d(seg_image)];
            image_to_show_rgb = zeros([size(image_to_show), 3]);
            image_to_show_rgb(:, :, 1) = image_to_show;
            image_to_show_rgb(:, :, 2) = binary_img_coordinates_2d;
            landmarks_round = round(landmarks);
            landmarks_img = zeros(size(seg_image));
            landmarks_img(landmarks_round(:, 2), landmarks_round(:, 1), landmarks_round(:, 3)) = 1;
            se = strel('square', 3);
            image_to_show_rgb(:, :, 3) = imdilate(reshape_2d(landmarks_img), se);
            figure, imshow(image_to_show_rgb, []);
            
        end
        seg_image_intensity = seg_image .* raw_image;
        seg_image_intensity = seg_image_intensity ./ sum(seg_image_intensity(:));
        feature = interp3(seg_image_intensity, Xq, Yq, Zq);
    case 'sphere'
        % convert to right-hand side, xyz -> yxz
        synapse_r= synapse;
        synapse_r(1) = synapse(2);
        synapse_r(2) = synapse(1);
        T = param.thickness;
        R = 8;
        [sphere_coordinates] = half_sphere_coordinates_extract(synapse_r, R);
        Xq = sphere_coordinates(:, 2);
        Yq = sphere_coordinates(:, 1);
        Zq = sphere_coordinates(:, 3);
        
        % visualize choosing region
        if false
            close all;
            round_coordinates = round(sphere_coordinates);
            binary_img_coordinates = zeros(size(seg_image));
            binary_img_coordinates(round_coordinates(:, 1), round_coordinates(:, 2), round_coordinates(:, 3)) = 1;
            
            binary_img_coordinates_2d = reshape_2d(binary_img_coordinates);
            % se = strel('diamond', 2);
            % binary_img_coordinates_2d = imdilate(binary_img_coordinates_2d, se);
            image_to_show = [];
            image_to_show = [image_to_show; reshape_2d(seg_image)];
            image_to_show_rgb = zeros([size(image_to_show), 3]);
            image_to_show_rgb(:, :, 1) = image_to_show;
            image_to_show_rgb(:, :, 2) = binary_img_coordinates_2d;
            landmarks_round = round(landmarks);
            landmarks_img = zeros(size(seg_image));
            landmarks_img(landmarks_round(:, 2), landmarks_round(:, 1), landmarks_round(:, 3)) = 1;
            se = strel('square', 3);
            % image_to_show_rgb(:, :, 3) = imdilate(reshape_2d(landmarks_img), se);
            figure, imshow(image_to_show_rgb, []);
            
        end
        seg_image_intensity = seg_image .* raw_image;
        seg_image_intensity = seg_image_intensity ./ sum(seg_image_intensity(:));
        feature = interp3(seg_image_intensity, Xq, Yq, Zq);        
              
end


end

function [cylinder_coordinates] = cylinder_coordinates_extract(origin, axis, T, R)

N = 10;
N_1 = 5;
N_2 = 5;
t = linspace(0, 2 * pi, N + 1);
t = t(1 : N);
% y: d1, x: d2
x = sin(t);
y = cos(t);

X = linspace(R, 1, N_1)' * x;
Y = linspace(R, 1, N_1)' * y;

X = repmat(X', [1, 1, N_2]);
Y = repmat(Y', [1, 1, N_2]);
Z = repmat(reshape(linspace(T + 1, 1, N_2), 1, 1, []), [N, N_1, 1]);

coordinates = [Y(:), X(:), Z(:), ones(numel(X), 1)];
orig_axis = [0, 0, 1];

translate_vec = origin - [0, 0, 0];
RM = compute_rotation_matrix_from_axises(orig_axis, axis);
TM = [RM, translate_vec'; 0, 0, 0, 1];
new_coordinates = TM * coordinates';
cylinder_coordinates = new_coordinates(1 : 3, :)';

% [TRI,v]= surf2patch(X,Y,Z,'triangle');
% 
% meshFV.faces = TRI;
% meshFV.vertices = v + 10;
% 
% rasterization_oversampling_scale = 4;
% segmentation_rasterization_function = @(given_mesh, window_size)rasterize_mesh(given_mesh, struct('oversampling_scale', rasterization_oversampling_scale, 'cropping', [ones(3, 1), window_size']));
% 
% I = segmentation_rasterization_function(meshFV, im_size);
% 
% uniq_x = unique(X(1, :));
% uniq_y = unique(Y(1, :));
% uniq_z = unique(Z(:, 1));

end

function [sphere_coordinates] = half_sphere_coordinates_extract(origin, R)

N = 5;
N_1 = 5;
N_2 = 5;
t = linspace(0, 2 * pi, N + 1);
t = t(1 : N);

phi = linspace(0, pi / 2, N_1 + 1);
phi = phi(2 : N_1 + 1);

x = bsxfun(@times, cos(t), sin(phi)');
z = bsxfun(@times, sin(t), sin(phi)');
y = bsxfun(@times, ones(1, numel(t)), cos(phi)');

R_n = linspace(1, R, N_2);

X = arrayfun(@(r)  r * x, R_n, 'UniformOutput', false);
X = cat(3, X{:});

Y = arrayfun(@(r)  r * y, R_n, 'UniformOutput', false);
Y = cat(3, Y{:});

Z = arrayfun(@(r)  r * z, R_n, 'UniformOutput', false);
Z = cat(3, Z{:});

coordinates = [Y(:), X(:), Z(:)];

sphere_coordinates = bsxfun(@plus, coordinates, origin);

end



function [rotation_matrix] = compute_rotation_matrix_from_axises(v1, v2)
% v1 : original axis, v2 : rotated axis.
% https://en.wikipedia.org/wiki/Rotation_matrix

v1 = v1(:);
v2 = v2(:);
v1 = v1 / norm(v1);
v2 = v2 / norm(v2);

% normal vector
n = cross(v1, v2);
% cos(theta) and sin(theta)
c = v1' * v2;
s = norm(n);

rotation_matrix = [c + n(1) ^ 2 * (1 - c), n(1) * n(2) * (1 - c) - n(3) * s, n(1) * n(3) * (1 - c) + n(2) * s; 
     n(2) * n(1) * (1 - c) + n(3) * s, c + n(2) ^ 2 * (1 - c), n(2) * n(3) * (1 - c) - n(1) * s; 
     n(3) * n(1) * (1 - c) - n(2) * s, n(3) * n(2) * (1 - c) + n(1) * s, c + n(3) ^ 2 * (1 - c)];

end







