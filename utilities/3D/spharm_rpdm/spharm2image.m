editfunction img = spharm2image(deg,fvec,meshtype,imageSize)

    if ~exist('meshtype','var')
        meshtype = [];
        meshtype.type = 'even';
    end
    mesh = [];
    [mesh.vertices,mesh.faces] = spharm2meshfigure(deg,fvec,meshtype);
    if ~exist('imageSize','var') 
        imageSize = ceil(max(mesh.vertices));
        imageSize = imageSize([2, 1, 3]);
    end
    img = surface_mesh_to_volume_image_conversion(mesh, imageSize);
    img = img > 0.5;
end
