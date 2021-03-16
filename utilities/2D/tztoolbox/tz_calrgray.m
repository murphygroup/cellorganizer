function img2=tz_calrgray(img,para,isinv)
%TZ_CALRGRAY Calibrate gray level of pixels according to positions.
%   IMG2 = TZ_CALRGRAY(IMG,PARA,ISINV) returns an image that is the
%   calibration of the [image] IMG. The calibartion is to adjust the image
%   intensity according to the postion relative to edge. The paramters of
%   calibration is PARA, which is a row vection with length 3. See
%   TZ_PROJBALL for more details.If ISINV is set to 1, the function will
%   try the inverse procedure. In IMG the pixels with intensity no greater
%   than 0 will be taken as background.
%   
%   See also

%   ??-???-2004 Initial write TINGZ
%   03-NOV-2004 Modified TINGZ
%       - add comments
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

img2=tz_edgedist(img);
distvec=tz_mat2vec(img2);
imgvec=tz_mat2vec(img);

distvec=(max(distvec)-distvec)/max(distvec);
egray=tz_projball(para,distvec);

eratio=max(egray)./egray(egray>0);
if isinv==1
    eratio=1./eratio;
end

imgvec(egray>0)=imgvec(egray>0).*eratio;

img2=reshape(imgvec,size(img));
