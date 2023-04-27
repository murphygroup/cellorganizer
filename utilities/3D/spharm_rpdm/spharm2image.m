function img = spharm2image(deg,fvec,options)
% edits
% 1/31/2023 R.F.Murphy correct default meshtype from "triangle" to "triangular"
% 3/23/2023 R.F.Murphy add option to control output image size
% 4/13/2023 R.F.Murphy change default debug option to false; fix oversampling_scale option
% 4/26/2023 R.F.Murphy flip and rotate image to match original orientation

default_options = [];
default_options.meshtype.type = 'triangular'; %1/31/2023
default_options.meshtype.nVertices = 4002;
default_options.oversampling_scale = 2; %4/13/2023
default_options.debug = false;
default_options.imagesize = [0, 0, 0];

if ~exist('options', 'var')
    options = default_options;
else
    options = process_options_structure(default_options, options);
end

mesh = [];
[mesh.vertices,mesh.faces] = spharm2meshfigure(deg,fvec,options.meshtype);
if max(options.imagesize)>0
    imageSize = options.imagesize;
    adjustLater = false;
else
    imageSize = ceil(max(mesh.vertices));
    imageSize = imageSize([2, 1, 3]);
    adjustLater = true;
end
img = surface_mesh_to_volume_image_conversion(mesh, imageSize, options);
imgf = flip(img,1); imgr = rot90(imgf,3); %4/26/2023
img = imgr > 0.5;
if adjustLater
    zmaxes = squeeze(max(max(img,[],1),[],2));
    nonemptyslices = find(zmaxes>0);
    img = img(:,:,min(nonemptyslices):max(nonemptyslices));
end
