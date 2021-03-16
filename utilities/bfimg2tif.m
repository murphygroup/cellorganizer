function answer = bfimg2tif( imgs, filename )

%xyzct
temp = zeros([size(imgs,1), size(imgs,2), 1, size(imgs,3), 1]);
for x=1:1:size(imgs,1)
    for y=1:1:size(imgs,2)
        for c=1:1:size(imgs,3)
            temp(x,y,1,c,1) = imgs(x,y,c);
        end
    end
end