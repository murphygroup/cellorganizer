function [ imgcompress, imsize, croprange ] = diffeo_compress_img( img )

    imsize = size(img);

    imsize = [imsize, ones(1,3 - length(imsize))];
    %Add extra zeros on the end incase the image is not 3D -grj

    %crop image to remove whitespace
    [cropimg, croprange] = cropImgND(img, 0); % use zero padding
    
    imgcompress = CompressLib.compress(cropimg);

end

