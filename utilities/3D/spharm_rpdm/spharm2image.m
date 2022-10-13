function img = spharm2image(deg,fvec,options)

default_options = [];
default_options.meshtype.type = 'triangle';
default_options.meshtype.nVertices = 4002;
default_options.options.oversampling_scale = 2;
default_options.debug = true;

if ~exist('options', 'var')
    options = default_options;
else
    options = process_options_structure(default_options, options);
end

mesh = [];
[mesh.vertices,mesh.faces] = spharm2meshfigure(deg,fvec,options.meshtype);
if ~exist('imageSize','var') 
    imageSize = ceil(max(mesh.vertices));
    imageSize = imageSize([2, 1, 3]);
end
img = surface_mesh_to_volume_image_conversion(mesh, imageSize, options);
img = img > 0.5;

end
