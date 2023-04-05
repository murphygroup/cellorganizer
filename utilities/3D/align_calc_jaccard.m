function jaccard_index = align_calc_jaccard(img1,img2)
% created 3/23/2023 R.F.Murphy

% requires that the two images are the same size
if(size(img1)~=size(img2))
    error('In jaccard_index: mismatched image sizes')
end

s1 = regionprops(img1,'centroid');
centroid1=s1.Centroid;
s2 = regionprops(img2,'centroid');
centroid2=s2.Centroid;
shift = round(centroid2-centroid1)
img1shifted = imtranslate(img1,shift);
centroid1shifted = regionprops(img1shifted,'centroid');
xorimg=xor(img1shifted,img2);
andimg=and(img1shifted,img2);
jaccard_index=sum(xorimg(:))/sum(andimg(:));