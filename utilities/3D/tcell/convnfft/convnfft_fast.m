function result = convnfft_faster( A, B )
% CONVNFFT_FASTER FFT-BASED N-dimensional convolution.
%   C = CONVNFFT(A, B) performs the N-dimensional convolution of
%   matrices A and B. Assumes that B is the kernel and returns a
%   matrix the same size as A.
%
% Class support for inputs A,B:
% float: double, single
%
% METHOD: CONVNFFT uses Fourier transform (FT) convolution theorem, i.e.
%         FT of the convolution is equal to the product of the FTs of the
%         input functions.
%         In 1-D, the complexity is O((na+nb)*log(na+nb)), where na/nb are
%         respectively the lengths of A and B.
%
% Usage recommendation:
%         In 1D, this function is faster than CONV for nA, nB > 1000.
%         In 2D, this function is faster than CONV2 for nA, nB > 20.
%         In 3D, this function is faster than CONVN for nA, nB > 5.
% 
% See also conv, conv2, convn.
% 
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original: 21-Jun-2009
%       23-Jun-2009: correct bug when ndims(A)<ndims(B)
%       02-Sep-2009: GPU/JACKET option
%       04-Sep-2009: options structure
%       16-Sep-2009: inplace product
%       11-Oct-2012: T. Buck removed GPU and Power2Flag options and the dims argument for now so the method can use fftn/ifftn.
%       5-Nov-2012: T. Buck removed all options except signal and kernel
%       2012-11-07 tebuck: renamed to convnfft_faster to differentiate with my computer's version of convnfft_fast.
%       2012-12-08 tebuck: allows odd-sized kernels (by padding the kernel by one towards the origin).

nd = max(ndims(A),ndims(B));
dims = 1:nd;

% This is non-periodic:
%A_freq = fftn(padarray(A, size(A) / 2));
%tebuck 11/8/2012
% whos
% Deal with odd-sized B (test this with something like convn([0, 0, 0, 1, 0, 0, 0], [1, 1, 1], 'same')):
B = padarray(B, mod([size(B), ones(1, 3 - ndims(B))], 2), 'pre');
A_freq = fftn(padarray(A, size(B) / 2));

% B_freq = fftn(circshift(padarray(B, size(B) / 2), size(B) - 1));
%tebuck 11/8/2012
%B_freq = fftn(circshift(padarray(B, size(A) / 2), size(A) - 1));
B_freq = fftn(circshift(padarray(B, size(A) / 2), (size(A) + size(B)) / 2 - 1));

result_freq = A_freq.*B_freq;
result = ifftn(result_freq);
% size_limits = [size(B) / 2; size(B) * 3 / 2 - 1];

%tebuck 11/8/2012
% size_limits = [size(A) / 2; size(A) * 3 / 2 - 1];
%size_limits = [size(B) / 2; size(B) * 3 / 2 - 1];
%result = result(size_limits(1, 1):size_limits(2, 1), size_limits(1, 2):size_limits(2, 2), size_limits(1, 3):size_limits(2, 3));

result_freq = A_freq.*B_freq;
result = ifftn(result_freq);
% size_limits = [size(B) / 2; size(B) * 3 / 2 - 1];
% size_limits = [size(A) / 2; size(A) * 3 / 2 - 1];
size_limits = [size(B) / 2; size(A) + size(B) / 2 - 1];
result = result(size_limits(1, 1):size_limits(2, 1), size_limits(1,2):size_limits(2, 2), size_limits(1, 3):size_limits(2, 3));
