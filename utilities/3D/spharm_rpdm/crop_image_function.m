function [I_out, start_ind] = crop_image_function(I_in, out_size)
% crop the image evenly along each edge and leave the center of the original image
% and also record the start indexs of the original image in the
% new image

if any(size(I_in) < out_size)
    error('pad size must be smaller or equal then the orignal image size')
end

if ~ismatrix(I_in) && ndims(I_in) ~= 3
    error('The dimension of image must be 2 or 3!')
end

img_size = size(I_in);

crop_size = floor((img_size - out_size) / 2);
start_ind = crop_size + 1;
end_ind = crop_size + out_size;

if ismatrix(I_in)
    I_out = I_in(start_ind(1) : end_ind(1), start_ind(2) : end_ind(2));
elseif ndims(I_in) == 3
    I_out = I_in(start_ind(1) : end_ind(1), start_ind(2) : end_ind(2), start_ind(3) : end_ind(3));
end
    
end