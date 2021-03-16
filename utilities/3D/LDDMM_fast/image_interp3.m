function c = image_interp3(a, b, method, extrapval, ~, ~)
% interp3(a, b, method, extrapval, pad_radius, window_mode)
% Finds values of a at locations b like Matlab's interp3.

if nargin < 3
    method = '*linear';
end
if nargin < 4
    extrapval = nan;
end

if iscell(b) && numel(b) == 3 && sum(abs(size(a) - size(b{1})))
    error('Images must be the same size')
end

if ismatrix(a)
    % xruan 09/04/2015 try c++ version interp2
    % vgg_interp2(varargin)
    %tic;
    if parallel.gpu.GPUDevice.isAvailable
        c = interp2(a, b{1}, b{2}, 'linear' );
    else
        c = vgg_interp2(a, b{1}, b{2}, method, extrapval);
    end
elseif ndims(a) == 3
    % xruan 10/31/2015 use mirt3D_mexinterp which is a 3D linear
    % interpolation method and is about 6-7 fold faster.
    % tic
    c = mirt3D_mexinterp(a, b{1}, b{2}, b{3});
    % toc
    % tic
    % c1 = interp3( a, b{1}, b{2}, b{3}, method, extrapval);
    % toc 
end
end
  
  
