function img2 = ml_improc(img,param)
%ML_IMPROC Process image.
%   IMG2 = ML_IMPROC(IMG,PARAM) returns an image that is the processed
%   version of the [image] IMG. PARAM is the parameters and has the
%   following fields:
%       'procfun' - processing function. This is a [general function].
%       'showmode' - show the result or not. 1 for show, 2 for both, 0 for
%           not.
%   
%   See also

%   14-Aug-2006 Initial write T. Zhao
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


if nargin < 2
    error('Exactly 2 arguments are required');
end

param = ml_initparam(param,struct('showmode',0));

img2 = ml_evalfun(img,param.procfun);

switch param.showmode
    case 1
        imshow(img2,[]);
    case 2
        subplot(1,2,1);
        imshow(img,[]);
        subplot(1,2,2);
        imshow(img2,[]);
end

