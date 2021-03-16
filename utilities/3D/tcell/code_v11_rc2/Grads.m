function gradmag = Grads( img1 )
if ndims(img1) == 3
    I = rgb2gray(img1);
else
    I = img1;
end
I=medfilt2(I); 

%Gauss filter
sigma=2; 
window=double(uint8(3*sigma)*2+2); 
H=fspecial('gaussian', window, sigma);
I=imfilter(I,H,'replicate');  

%Gradient segmentation
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

gradmag(gradmag>20)=255;
gradmag(gradmag<=20)=0;

end

