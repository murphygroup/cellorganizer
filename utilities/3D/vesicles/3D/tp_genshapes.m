function [nucimg,cellimg] = tp_genshapes(index)
% Generate nuclear and cell shapes

load('./inter_results/nuc_shape_spline/nuc_shape_model')
disp('Generating nuclei shape ...')
instance = tp_genspsurf(model);
tp_gennucshape_result = tp_gennucshape(instance);
nucimg = tp_gennucshape_result.nucimg;
nucsurf = tp_gennucshape_result.nucsurf;
nuclei.nucimgsize = size(nucimg);
nuclei.nucsurf = nucsurf;

load('./inter_results/cell_shape_eigen/cell_shape_model')
disp('Generating cell shape ...')
cellimg = tp_gencellshape(model,nuclei);

disp('Resizing cell shape ...')
box = tp_imbox(cellimg);
nucimg = nucimg(box(1):box(2),box(3):box(4),:);

cellimg = cellimg(box(1):box(2),box(3):box(4),:);
cellimg = tp_stretch3d(cellimg,nuclei.nucimgsize(3));

disp('Saving ...')
blank = nucimg(:,:,1);
cellimg = cat(3,blank,cellimg,blank);
nucimg = cat(3,blank,nucimg,blank);

save(['./inter_results/image_synshapes/shape'...
    int2str(index) '.mat'], 'nucimg', 'cellimg');
