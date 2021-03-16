function shape = ml_mxs2mxp(medaxis,width,param)
%ML_MXS2MXP Convert medial axis reprenstation into spline representation.
%   SHAPE = ML_MXS2MXP(MEDAXIS,WIDTH) returns a structure that is the 'mxp'
%   shape description. MEDAXIS is the medial axis and WIDTH is the width.
%   
%   SHAPE = ML_MXS2MXP(MEDAXIS,WIDTH,PARAM) lets the user specify
%   parameters for splines. See TZ_MEDAXSPFEAT for more details.
%
%   See also TZ_MEDAXSPFEAT

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

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

if nargin < 2
    error('Exactly 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

shape.format = 'mxp';

feat = ml_medaxspfeat(medaxis,width,param);

shape.length = feat{1};
shape.spmedaxis = ml_feat2sp(feat{2},feat{3});
shape.spwidth = ml_feat2sp(feat{4},feat{5});

