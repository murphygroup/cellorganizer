function [ img ] = gauss2img( centers, covars, objinten, imsize, thresh )
%centers - nx3 matrix of gaussian object centesr
%covars - 3x3xn matrix of gaussian object covariance matrices
%objinten - nx1 vector of object intensities
%imgsize - 1x3 vector dimensions of image
%thresh - optional - cutoff gaussian objects at some threshold

if ~exist('thresh', 'var') | isempty(thresh)
    thresh = false;
end

img = zeros(imsize);

for i = 1:size(centers,1)
    try
        obj = ml_gaussimg(covars(:,:,i)) * objinten(i);

        if thresh
            prct = prctile(obj(:), thresh);
            obj = double(obj >= prct);
            obj = placeObject(obj, imsize, centers(i,:)) >= 1;
            img(obj) = 1;
        else

            img = img + placeObject(obj, imsize, centers(i,:));
        end

    catch
        disp(['Could not place object ' num2str(i)]);
    end
end


end

