function [kernel] = Greedy3D_kernel_resampled(s, a, g, use_second_version)
  % Internal function for Greedy3D_lambda_pre_compressed.
  % 
  % This function returns the low pass filter, i.e., the squared inverse of the differential operator paper from Rohde et al. 2008 ((alpha * del^2 + gamma)^-2) in the spatial domain.
  %
  % 2011-06-02 17:25 tebuck: double continuous Laplacian using 
  % 2011-06-02 18:02 tebuck: Testing by hand indicates that the
  % radius of the discrete kernel is about 64 pixels. This function
  % will return resampled versions of the kernel with the scale of
  % one set to 128 x 128 x [input z size].
  % 2011-06-02 20:16 tebuck: Scale of one set to 64 x 64 x [input z
  % size].
  % 2011-07-07 10:17 tebuck: It seems to me that the first version
  % improperly moves the center of the kernel toward the actual
  % center of the image (between the center four pixels with even
  % sizes), and the second version doesn't interpolate between the
  % edges of the image, which will end up being adjacent to each
  % other... Maybe I should switch to having odd-sized kernels. I
  % think convnfft can handle them.
  % 2011-08-06 17:15 tebuck: Previous changes led to the second
  % version using a kernel that is 64 cubed before scaling. Current
  % change is that rescale_kernel_by_sum controls whether the kernel
  % is normalized to the same sum it had prior to scaling. Hopefully
  % this will prevent a smaller kernel from producing large
  % velocities. It's off by default, so the behavior should remain
  % the same unless the option is used.
  % It just occurred to me that the renormalization would give any
  % preference to larger dimensions because of the way imresize
  % works, so this option is unnecessary.

if ~exist('use_second_version', 'var')
  use_second_version = false;
end

s(numel(s) + 1:3) = 1; 

M = 64; 
N = 64; 
P = max(s(3), 1); 
if use_second_version
  P = 64;
end
M2 = s(1); 
N2 = s(2); 
P2 = max(s(3), 1); 
kkx = (2*pi/N)*[0:(N/2-1) (-N/2):(-1)];
kky = (2*pi/M)*[0:(M/2-1) (-M/2):(-1)];
kkz = (2*pi/P)*[0:(P/2-1) (-P/2):(-1)];
% Make a 2D kernel if the Z size is 1 or absent:
if P2 == 1
  kkz = zeros(size(kkz));
end

[KX,KY,KZ] = meshgrid(kkx,kky,kkz); 
delsq = 2*a*(...
  (1-cos(KX)) + ...
  (1-cos(KY)) + ...
  (1-cos(KZ)) ...
  ) + g;
delsq = delsq * ((M/M2));
delsq = -(delsq.^2);
delsq(1,1,1) = 1;
delsq = delsq .* exp(...
  i * KX * N * .5 + ...
  i * KY * M * .5 + ...
  i * KZ * P * .5) ;
kernel = ifftn(1 ./ delsq);
kernel = real(kernel);
kernel_sum = abs(sum(kernel(:)));
kernel = imresize(kernel, [M2, N2], 'bilinear');
if use_second_version
  % 3D rescale hack:
  kernel = permute(kernel, [1, 3, 2]);
  kernel = imresize(kernel, [M2, P2], 'bilinear');
  kernel = ipermute(kernel, [1, 3, 2]);
end
kernel = kernel .* kernel_sum ./ abs(sum(kernel(:)));
