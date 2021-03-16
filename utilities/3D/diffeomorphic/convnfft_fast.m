function result = convnfft_faster(A, B, padding_method, corr)
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
%       2013-01-04 tebuck: allows odd-sized A as well as B by the same mechanism. Tests:
%       (convn([0, 0, 0, 1, 0, 0], [1, 1], 'same') - convnfft_fast([0, 0, 0, 1, 0, 0], [1, 1]))
%       (convn([0, 0, 0, 1, 0, 0], [1, 1, 1], 'same') - convnfft_fast([0, 0, 0, 1, 0, 0], [1, 1, 1]))
%       (convn([0, 0, 0, 1, 0, 0, 0], [1, 1, 1], 'same') - convnfft_fast([0, 0, 0, 1, 0, 0, 0], [1, 1, 1]))
%       (convn([0, 0, 0, 1, 0, 0, 0], [1, 1], 'same') - convnfft_fast([0, 0, 0, 1, 0, 0, 0], [1, 1]))
%       2013-03-13 tebuck: allows selecting the padding method for the first input.


  if ~exist('padding_method', 'var') || isempty(padding_method)
    % Hack to get the default zero-fill behavior:
    padding_method = 0;
  end
  
  if ~exist('corr', 'var') || isempty(corr)
      corr = false;
  end

  A_size = size(A);
  B_size = size(B);

  nd = max(ndims(A),ndims(B));


  % nd, dims, dims_logical

  dims = 1:nd;
  % Pad to make this non-periodic:
  % Deal with odd-sized A:
  A_size_odd = mod([size(A), ones(1, 3 - ndims(A))], 2);
  % A_size_odd = mod([size(A), ones(1, nd - ndims(A))], 2);
  A = padarray(A, A_size_odd, padding_method, 'pre');
  % Deal with odd-sized B:
  B = padarray(B, mod([size(B), ones(1, 3 - ndims(B))], 2), 'pre');
  % B = padarray(B, mod([size(B), ones(1, nd - ndims(B))], 2), 'pre');

  % Pad A and B prevent interactions of opposite sides of images due to the transform's periodicicity:
  A_freq = fftn(padarray(A, size(B) / 2, padding_method));
  B_freq = fftn(circshift(padarray(B, size(A) / 2), (size(A) + size(B)) / 2 - 1));
  
  if ~corr
    result_freq = A_freq.*B_freq;
  else
    result_freq = A_freq.*conj(B_freq);
  end
  % tic
  result = ifftn(result_freq);
  % toc

  size_limits = [size(B) / 2; size(A) + size(B) / 2 - 1];
  % Deal with odd-sized A:
  % keyboard
  size_limits(1, :) = size_limits(1, :) + A_size_odd;
  result = result(size_limits(1, 1):size_limits(2, 1), size_limits(1,2):size_limits(2, 2), size_limits(1, 3):size_limits(2, 3));

