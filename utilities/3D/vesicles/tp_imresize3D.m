function img2 = tp_imresize3D(img,zres,method)

%TP_IMRESIZE3D resize 3D iamges to uniformize the resolutions

zdim = size(img,3);

switch method
    case 'downsampling'
        for z = 1:zdim
            img2(:,:,z) = imresize(img(:,:,z),1/zres);
        end
    case 'upsampling'
        img2 = interp3(img);
end