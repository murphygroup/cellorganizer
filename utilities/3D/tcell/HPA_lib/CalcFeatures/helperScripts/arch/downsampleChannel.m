function img = downsampleChannel( origImg, dsf, img)

s = size(origImg);
if isempty(img)
    img = imresize(origImg,[s(1) s(2)]/dsf);
end

return

