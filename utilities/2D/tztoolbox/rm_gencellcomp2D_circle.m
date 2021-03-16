function [nucimg,cellimg] = rm_gencellcomp2D_circle(model,param)

imageSize=param.imageSize(1);
xymax = 1.25;
inc = xymax/((imageSize-1)/2);
[meshx,meshy]=meshgrid(-xymax:inc:xymax);
meshx=reshape(meshx,size(meshx,1)^2,1);
meshy=reshape(meshy,size(meshy,1)^2,1);
radii = (meshx.^2+meshy.^2);
eps=0.003;
idxcell = abs(1-radii)<eps;
idxnuc = abs(0.25-radii)<eps;

nucimg = zeros(imageSize);
cellimg = nucimg;

nucimg(idxnuc) = 1;
cellimg(idxcell) = 1;



