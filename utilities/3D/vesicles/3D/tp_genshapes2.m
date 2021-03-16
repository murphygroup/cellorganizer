function [nucimg,cellimg] = tp_genshapes2(model)

disp('Generating nuclear shape ...')
instance = tp_genspsurf(model.nuclearShapeModel);
tp_gennucshape_result = tp_gennucshape(instance);
nucimg = tp_gennucshape_result.nucimg;
nucsurf = tp_gennucshape_result.nucsurf;
nuclei.nucimgsize = size(nucimg);
nuclei.nucsurf = nucsurf;

disp('Generating cell shape ...')
cellimg = tp_gencellshape(model.cellShapeModel,nuclei);

disp('Resizing cell shape ...')
box = tp_imbox(cellimg);
nucimg = nucimg(box(1):box(2),box(3):box(4),:);
cellimg = cellimg(box(1):box(2),box(3):box(4),:);
cellimg = tp_stretch3d(cellimg,nuclei.nucimgsize(3));

disp('Saving ...')
blank = nucimg(:,:,1);
cellimg = cat(3,blank,cellimg,blank);
nucimg = cat(3,blank,nucimg,blank);
