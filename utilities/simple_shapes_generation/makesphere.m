function S = makesphere(imgsize, r);
[xx, yy, zz] = meshgrid(1:imgsize, 1:imgsize, 1:imgsize);
if mod(imgsize, 2) == 0
    c = imgsize/2;
else
    c = ceil(imgsize/2);
end
S = sqrt((xx-c).^2+(yy-c).^2+(zz-c).^2)<=(r);
end%makesphere
