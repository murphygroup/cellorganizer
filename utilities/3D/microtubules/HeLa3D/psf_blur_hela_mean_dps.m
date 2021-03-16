function G_psf = psf_blur_hela_mean_dps(G,psf)
%PSF_BLUR_HELA_MEAN_DPS

% Author: Jieyue Li and Devin Sullivan
% Edited: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
%
% Copyright (C) 2011-2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

if nargin ~= 2
   warning('CellOrganizer: Wrong number of input arguments');
end

G_psf = [];

if isempty(psf)
  warning('CellOrganizer: Input argument psf is empty. Returning original image');
  G_psf = G;
  return
else
   try
      G_psf = imfilter(G,psf,'conv','same');
   catch
      warning('CellOrganizer: Failed to convolve image with given point spread function. Returning original image');
      G_psf = G;
      return
   end
end
