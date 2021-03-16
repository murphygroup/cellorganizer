function I=imgaussian(I,sigma,siz)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an 
% Gaussian filter. This function uses for filtering IMFILTER or if 
% compiled the fast  mex code imgaussian.c . Instead of using a 
% multidimensional gaussian kernel, it uses the fact that a Gaussian 
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: The 1D, 2D greyscale/color, or 3D input image with 
%           data type Single or Double
%   SIGMA: The sigma used for the Gaussian kernel
%   SIZE: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   J: The gaussian filtered image
%
% note, compile the code with: mex imgaussian.c -v
%
% example,
%   I = im2double(imread('peppers.png'));
%   figure, imshow(imgaussian(I,10));
% 
% Function is written by D.Kroon University of Twente (September 2009)
% 
% 2012-12-07 tebuck: using convnfft_fast to reduce running time for larger sigma.

% filter_function = @(x, y)imfilter(x, y, 'same' ,'replicate');
filter_function = @(x, y)convnfft_fast(x, y);

if(~exist('siz','var')), siz=sigma*6; end

if(sigma>0)
    % Make 1D Gaussian kernel
    x=-ceil(siz/2):ceil(siz/2);
    H = exp(-(x.^2/(2*sigma^2)));
    H = H/sum(H(:));

    % Filter each dimension with the 1D Gaussian kernels\
    if(ndims(I)==1)
        I=filter_function(I,H);
    elseif(ndims(I)==2)
        Hx=reshape(H,[length(H) 1]);
        Hy=reshape(H,[1 length(H)]);
        I=filter_function(filter_function(I,Hx),Hy);
    elseif(ndims(I)==3)
        if(size(I,3)<4) % Detect if 3D or color image
            Hx=reshape(H,[length(H) 1]);
            Hy=reshape(H,[1 length(H)]);
            for k=1:size(I,3)
                I(:,:,k)=filter_function(filter_function(I(:,:,k),Hx),Hy);
            end
        else
            Hx=reshape(H,[length(H) 1 1]);
            Hy=reshape(H,[1 length(H) 1]);
            Hz=reshape(H,[1 1 length(H)]);
            I=filter_function(filter_function(filter_function(I,Hx),Hy),Hz);
        end
    else
        error('imgaussian:input','unsupported input dimension');
    end
end