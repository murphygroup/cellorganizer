function raw2proc_hela(cellnum)

protim3 = ml_loadimage(['raw/cell' num2str(cellnum) '/prot/'],'tif');
dnaim3 = ml_loadimage(['raw/cell' num2str(cellnum) '/dna/'],'tif');
cellim3 = ml_loadimage(['raw/cell' num2str(cellnum) '/cell/'],'tif');
mask = ml_loadimage(['raw/cell' num2str(cellnum) '/crop/'],'tif');

protim3 = tz_maskimg_3d(protim3,mask);
dnaim3 = tz_maskimg_3d(dnaim3,mask);
cellim3 = tz_maskimg_3d(cellim3,mask);

% Downsize the images
protim3 = ml_downsize(protim3,[4 4 1],'average');
dnaim3 = ml_downsize(dnaim3,[4 4 1],'average');
cellim3 = ml_downsize(cellim3,[4 4 1],'average');

% Resize the image so that resolutions match in X,Y,Z
mkdir(['proc/cell_' num2str(cellnum)]);
save(['proc/cell_' num2str(cellnum) '/cell' num2str(cellnum) '.mat']);
