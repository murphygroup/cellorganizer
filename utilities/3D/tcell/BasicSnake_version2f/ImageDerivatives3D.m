function J=ImageDerivatives3D(I,sigma,type, voxel_scales)
% Gaussian based image derivatives
%
%  J=ImageDerivatives3D(I,sigma,type)
%
% inputs,
%   I : The 3D image
%   sigma : Gaussian Sigma
%   type : 'x', 'y', 'z', 'xx', 'yy', 'zz', 'xy', 'xz', 'yz'
%
% outputs,
%   J : The image derivative
%
% Function is written by D.Kroon University of Twente (July 2010)
% 
% 2012-12-06 tebuck: using convnfft_fast to reduce running time for larger sigma.
% 2013-01-20 tebuck: added voxel_scales to account for non-uniform voxel dimensions (either weakening or strengthening by using the reciprocal).

% Make derivatives kernels
[x,y,z]=ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));

switch(type)
    case 'x'
        DGauss=-(x./((2*pi)^(3/2)*sigma^5)).*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2)) * voxel_scales(1);
    case 'y'
        DGauss=-(y./((2*pi)^(3/2)*sigma^5)).*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2)) * voxel_scales(2);
    case 'z'
        DGauss=-(z./((2*pi)^(3/2)*sigma^5)).*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2)) * voxel_scales(3);
    case 'xx'
        DGauss = 1/((2*pi)^(3/2)*sigma^5) * (x.^2/sigma^2 - 1) .* exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2)) * voxel_scales(1)^2;
    case 'yy'
        DGauss = 1/((2*pi)^(3/2)*sigma^5) * (y.^2/sigma^2 - 1) .* exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2)) * voxel_scales(2)^2;
    case 'zz'
        DGauss = 1/((2*pi)^(3/2)*sigma^5) * (z.^2/sigma^2 - 1) .* exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2)) * voxel_scales(3)^2;
    case {'xy','yx'}
        DGauss = 1/((2*pi)^(3/2)*sigma^7) * (x .* y)           .* exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2)) * prod(voxel_scales([1, 2]));
    case {'xz','zx'}
        DGauss = 1/((2*pi)^(3/2)*sigma^7) * (x .* z)           .* exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2)) * prod(voxel_scales([1, 3]));
    case {'yz','zy'}
        DGauss = 1/((2*pi)^(3/2)*sigma^7) * (y .* z)           .* exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2)) * prod(voxel_scales([2, 3]));
end

% % tic
% K=SeparateKernel(DGauss);
% % toc
% J = imfilter(I,K{1},'conv','symmetric');
% J = imfilter(J,K{2},'conv','symmetric');
% J = imfilter(J,K{3},'conv','symmetric');
% tic
% keyboard
J = convnfft_fast(I, DGauss);
% toc
