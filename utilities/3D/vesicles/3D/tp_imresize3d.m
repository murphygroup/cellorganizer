function img_resized = tp_imresize3d(img,newSize,method,max_memory)

% Resize a 3D stack image
%   img - input 3D image
%   newsize - the target image size
%   method - the interpolation method for resizing

if nargin < 4
    max_memory = 4e9;
end
img_memory = getfield(whos('img'), 'bytes');
img_sizeof = img_memory / numel(img);
img_resized_memory = prod(newSize) * img_sizeof;

% Original image coordinate
[Ydim, Xdim, Zdim] = size(img);
Yhalf = (Ydim-1)/2;
Xhalf = (Xdim-1)/2;
Zhalf = (Zdim-1)/2;
Xvalues = -Xhalf:Xhalf;
Yvalues = -Yhalf:Yhalf;
Zvalues = -Zhalf:Zhalf;

% New image coordinate
Yhalf2 = (newSize(1)-1)/2;
Xhalf2 = (newSize(2)-1)/2;
Zhalf2 = (newSize(3)-1)/2;
R = size(img)./newSize;     % The resizing ratio
Xvalues2 = (-Xhalf2:Xhalf2)*R(2);
Yvalues2 = (-Yhalf2:Yhalf2)*R(1);
Zvalues2 = (-Zhalf2:Zhalf2)*R(3);

% Interpolation
% img, X, Y, Z, passed img, X, Y, Z, guess two temporary copies, Xi, Yi, Zi, passed Xi, Yi, Zi, result img_resized, guess two temporary copies
required_memory = 10 * img_memory + 9 * img_resized_memory;
if required_memory <= max_memory
    [X,Y,Z] = meshgrid(-Xhalf:Xhalf,-Yhalf:Yhalf,-Zhalf:Zhalf);
    [Xi,Yi,Zi] = meshgrid(Xvalues2,Yvalues2,Zvalues2);
    img_resized = interp3(X, Y, Z, img, Xi, Yi, Zi, method);
else
    % Multiple interpolations to be more memory efficient
    warning('Assumes interpolation method is separable to be more memory efficient')
    [X,Y] = meshgrid(Xvalues,Yvalues);
    [Xi,Yi] = meshgrid(Xvalues2,Yvalues2);
    img_resized = zeros([newSize(1), newSize(2), Zdim], class(img));
    for i = 1:Zdim
        img_resized(:, :, i) = interp2(X, Y, img(:, :, i), Xi, Yi, method);
    end
    dimensions_permutation = [3, 2, 1];
    img_resized = permute(img_resized, dimensions_permutation);
    img_resized2 = zeros(newSize(dimensions_permutation), class(img));
    [X2,Z] = meshgrid(1:newSize(2),Zvalues);
    [Xi2,Zi] = meshgrid(1:newSize(2),Zvalues2);
    for i = 1:newSize(1)
        img_resized2(:, :, i) = interp2(X2, Z, img_resized(:, :, i), Xi2, Zi, method);
    end
    img_resized = ipermute(img_resized2, dimensions_permutation);
end

