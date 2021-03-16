function segmentedregions(regions,prot,ifiles,inimages)

bwl = bwlabel(regions==max(regions(:)),4);
s = size(regions);
ma = max(bwl(:));
%figure
mkdir('./segmented_images',ifiles(inimages,1).name);  %make directory to store the segmented images
cd(strcat('./segmented_images/',ifiles(inimages,1).name));
for i=1:ma
    [r c] = find(bwl==i);
    minr = min(r);  maxr = max(r);  minc = min(c);  maxc = max(c);
    newIm = uint8(zeros(maxr-minr+1,maxc-minc+1));
    idx = (c-1)*s(1)+r;
    idx2 = (c-min(c))*(max(r)-minr+1)+r-min(r)+1;
    newIm(idx2) = prot(idx);
    %imshow(newIm);
    %pause
    
    imwrite(newIm,strcat(num2str(i),'.png'));
end
