function [img,imgfiles]=ml_loadimage(imgdir,ext,loadthre)
%ML_LOADIMAGE Read images under one directory into a 3D matrix.
%   IMG = ML_LOADIMAGE(IMGDIR,EXTENSION) returns a 3D matrix IMG by reading
%   images with extension EXTENSION from the directory IMGDIR. The image files
%   are sorted by their number labels before reading.
%
%   IMG = ML_LOADIMAGE(MGDIR,EXTENSION,LOADTHRESHOLD) set all values greater
%   than LOADTHRESHOLD to 0 in the 3D matrix. 
%   ML_LOADIMAGE(MGDIR,EXTENSION,[]) is the same as 
%   ML_LOADIMAGE(IMGDIR,EXTENSION).
%
%   [IMG,IMGFILES] = ML_LOADIMAGE(...) also returns sorted file names in
%   a cell array IMGFILES.

% Copyright (C) 2006  Murphy Lab
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

%   ??-???-???? Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add comments
%   15-MAR-2004 Modified TINGZ
%       - get strictly increasing file numbers
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin<2
    error('2 or 3 arguments are required.');
end

if nargin<3
    loadthre=[];
end

imgfiles=ml_dir([imgdir '/*.' ext]);

if(isempty(imgfiles))
    img=[];
    imgfiles={};
    return;
end

pos=1;

for i=1:length(imgfiles)
    filenum(i)=ml_getfilenum(imgfiles{i});
end
[sorted,num]=sort(filenum);

if length(sorted)>0
    while any(sorted(1:end-1)-sorted(2:end)==0) & all(num>0)
        pos=pos+1;
        for i=1:length(imgfiles)
            filenum(i)=ml_getfilenum(imgfiles{i},pos);
        end
        if any(filenum<0)
            break
        else
            [sorted,num]=sort(filenum);
        end
    end 
end

imgfiles=imgfiles(num);

for i=length(imgfiles):-1:1
    img(:,:,i)=ml_readimage([imgdir '/' imgfiles{i}]);
end

if ~isempty(loadthre)
    img(find(img > loadthre)) = 0;
end