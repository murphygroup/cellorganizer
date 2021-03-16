function seeds = bgseeds( img, MINNUCLEUSDIAMETER, IMAGEPIXELSIZE)

minarea = (MINNUCLEUSDIAMETER/IMAGEPIXELSIZE/2)^2*pi;

seeds = imerode(img==min(img(:)),strel('disk',2));
bw = seeds>0;
bwl = bwlabel(bw,4);
props = regionprops(bwl,'Area');
stats = zeros(size(props));
for i=1:length(stats)
    stats(i) = props(i).Area;
end
idx = find(stats<minarea);
bwmask = ismember(bwl,idx);
seeds(bwmask==1) = 0;

seeds = uint8(seeds);
