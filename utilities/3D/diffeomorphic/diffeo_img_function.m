function [img] = diffeo_img_function(filenum, compressed_imgs, maxsize, imsizes, imcrops)

% Taraz Buck and Greg Johnson
%
% Copyright (C) 2012-2015 Murphy Lab
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu

if ~isempty(compressed_imgs(filenum))
    %icaoberg 2015/12/16 just prints the image index
    disp( ['Decompressing image ' num2str(filenum)] );
    img = CompressLib.decompress(compressed_imgs{filenum});
    
    img = padarray(img, [imcrops(filenum,1)-1, imcrops(filenum,3)-1], 'pre');
    img = padarray(img, [imsizes(filenum,1) - imcrops(filenum,2), imsizes(filenum,2) - imcrops(filenum,4)], 'post');
    
else
    img = [];
end

if ~isempty(img)
    imsize = size(img);
    imsize = [imsize, ones(1,3 - length(imsize))];
    
    pad = (maxsize - imsize) ./ 2;
    
    img = padarray(img, floor(pad), 0, 'both');
    
    img = padarray(img, double(~(floor(pad) == pad)), 0, 'pre');
    img = double(img);
end
end