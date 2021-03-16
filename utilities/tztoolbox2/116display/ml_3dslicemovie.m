function mov = ml_3dslicemovie(img,param)
%ML_3DSLICEMOVIE Show slices of a 3D image in a movie.
%   MOV = ML_3DSLICEMOVIE(IMG) returns a matlab movie that is for the slices
%   of the 3D image IMG.
%   
%   MOV = ML_3DSLICEMOVIE(IMG,PARAM)
%   
%   See also

%   21-Oct-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('mode','3D','mask',[],'sizeratio',[]));

if ~isempty(param.mask)
    img = tz_maskimg_3d(img,param.mask);
end

switch param.mode
    case 'ID',
        no_of_slices=size(img,3);
        
        for slice_no = 1:no_of_slices
            slice = img(:,:,slice_no);
            if isempty(param.sizeratio)
                imshow(slice,double([min(slice(:)),max(slice(:))]));
            else
                imshow(slice,double([min(slice(:)),max(slice(:))]), ...
                       'InitialMagnification',param.sizeratio);
            end
            
            title(slice_no)
            mov(slice_no) = getframe;
        end
                
    case '3D',
        minpixel=min(img(:))
        maxpixel=max(img(:))
        no_of_slices=size(img,3);
        %[x,y,z]=find(img==maxpixel);
        for slice_no = 1:no_of_slices
            slice = img(:,:,slice_no);
            %[x,y]=find(slice==maxpixel);
            if isempty(param.sizeratio)
                imshow(slice,double([minpixel,maxpixel]));
            else
                imshow(slice,double([min(slice(:)),max(slice(:))]), ...
                       'InitialMagnification',param.sizeratio);
            end
            
            title(slice_no)
            mov(slice_no) = getframe;
        end
    otherwise,
        warning('Invalid mode');
end
