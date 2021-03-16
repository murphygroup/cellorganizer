function loc = MarkerControlled_Watershed2(HandImage)  
if size(HandImage,3)>1
    I = rgb2gray(HandImage);  
else
    I  =  HandImage;  
end
I=double(I);
I=I.^1.5;
I=floor(256*(I-min(min(I)))/(max(max(I))-min(min(I))));
se = strel('disk',5);
Ie = imerode(uint8(I), se); 
Iobr = imreconstruct(Ie, uint8(I));%first corrosion and then rebuild  
se3 = strel('disk',4); 
Iobrd = imdilate(Iobr, se3);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));  
Iobrcbr = imcomplement(Iobrcbr);
fgm = imregionalmax(Iobrcbr);  
I2 = I;  
I2(fgm) = 255;% fgm of foreground area (The pixel value is 1) mark to original image
se2 = strel(ones(4,4));  
fgm2 = imclose(fgm, se2);  
fgm3 = imerode(fgm2, se2);  
fgm4 = bwareaopen(fgm3, 10);  
I3 = I;  
I3(fgm4) = 255;  
bw=im2bw(Iobrcbr,graythresh(Iobrcbr));  
D = bwdist(bw);  
DL = watershed(D);  
bgm = DL == 0;%In the watershed transformation result L, the same area is represented by the same number, and the interregional boundary is identified by 0
%Calculate the longitudinal gradient
hy = fspecial('sobel');
hx = hy'; 
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');  
gradmag = sqrt(Ix.^2 + Iy.^2);%Calculate the gradient magnitude 
gradmag2 = imimposemin(gradmag, bgm | fgm4);  
L = watershed(gradmag2);  
I4 = I;  
I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;    

num=unique(L);
cnt=1;
for ll=1:length(num)
    [m,n]=find(L==num(ll));
    tmp=sum(L==num(ll));
    if tmp<=40
        loc{cnt,1}=[sum(n)/length(n),sum(m)/length(m)];
        % plot(sum(n)/length(n),sum(m)/length(m),'r*');
        cnt=cnt+1;
    end
end