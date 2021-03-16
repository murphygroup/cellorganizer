function OutImg=imadd_img(Img)
img_num = length(Img);
if size(Img{1,1},3)==3
    img=rgb2gray(Img{1,1});
else
    img=Img{1,1};
end
img_add=zeros(size(Img{1,1}));
for j = 1:img_num %Read images in turn
    Y_k =   Img{j,1};
    q2=size(Y_k,3);
    if q2>1
        Y_k=rgb2gray(Y_k);
    end
    Y_k=double(Y_k);
    img_add = imadd(img_add,Y_k) ;
end
OutImg = Normalize0_255(img_add) ;
OutImg=uint8(OutImg);