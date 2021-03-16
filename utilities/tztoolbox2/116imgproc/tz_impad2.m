function img3 = tz_impad2(img2,img1)
%TZ_IMPAD2 Pad or cut image to the same size of the other image
%   IMG3 = TZ_IMPAD2(IMG1,IMG2) simply pad or cut IMG1 to the size of IMG2.
 
%   ??-Oct-2005 Initial write T. Zhao

img3 = img2;

if size(img2,1) < size(img1,1)
    img3 = [img3;zeros(size(img1,1)-size(img3,1),size(img3,2))];
end

if size(img2,2) < size(img1,2)
    img3 = [img3,zeros(size(img3,1),size(img1,2)-size(img3,2))];
end

img3 = img3(1:size(img1,1),1:size(img1,2));