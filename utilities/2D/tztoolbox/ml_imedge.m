function img2=ml_imedge(img,cropimg,option)
%ML_IMEDGE Extract edges for objects in the image.
%   IMG2 = ML_IMEDGE(IMG,CROPIMG,OPTION) returns an image that contains the
%   edges of objects in the image IMAGE, which will be masked first by
%   CROPIMG before edge detection. If CROPIMG is empty([]), there will be
%   no masking. OPTION specifies how to do edge detection:
%       'bw' - use BWPERIM directly.
%       'ce' - cell edge
%       'nu' - nuclear edge
%   This function is especially designed for extracting edges of cells and 
%   nucleus.
%
%   See also

%   27-SEP-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
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

if nargin < 3
    error('Exactly 3 arguments are required')
end

img=double(img);

switch(option)
case 'bw'
    if ~isempty(cropimg)
        img(cropimg==0)=0;
    end
    img2=bwperim(img);
case 'ce'
    common = ml_imgcommonpixel2D(img);
    img2=bwmorph(img>common,'erode');
    if ~isempty(cropimg)
        img2(cropimg==0)=0;
    end
    img2 = bwperim(imfill(img2>0,'hole'));
%     medimg=medfilt2(img,[5 5]);
%     img2=bwperim(imfill(medimg>0,'hole'));
case 'nu'
    [img,maskimage]=ml_preprocess(double(img),cropimg,'ml','yesbgsub');
    img2=bwperim(imfill(img>0,'hole'));
end
    
