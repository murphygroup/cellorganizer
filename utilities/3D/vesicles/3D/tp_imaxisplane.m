function [medplane,height,mask,rbase,cbase] = tp_imaxisplane(img)
%ML_IMAXISPLANE Extract and process medial axis plane from an 3d image.
%   IMGAXIS = ML_IMAXIS(IMG) returns an image is superimposition of IMG and
%   its mecial axis. This function will take all pixels above intensity 0 
%   as objects.

zindex = repmat(shiftdim(1:size(img,3),-1),[size(img,1),size(img,2)]);
zindex(img==0) = NaN;
topsurf = min(zindex,[],3);
bottomsurf = max(zindex,[],3);
medplane = (topsurf + bottomsurf) / 2;
height = (bottomsurf - topsurf) / 2;

mask = medplane > 0;
rbase = find(any(mask,2));
cbase = find(any(mask,1));

medplane = medplane(rbase,cbase);
height = height(rbase,cbase);
mask = medplane > 0;

bwimg = imfill(img,'holes');
fluobylayer = squeeze(sum(sum(bwimg)));
[m,medind] = max(fluobylayer);

medplane(isnan(medplane)) = medind;
height(isnan(height)) = 0;