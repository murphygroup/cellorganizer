function shape = tz_mxs2mxp(medaxis,width,param)
%TZ_MXS2MXP Obsolete. See ML_MXS2MXP.
%   SHAPE = TZ_MXS2MXP(MEDAXIS,WIDTH) returns a structure that is the 'mxp'
%   shape description. MEDAXIS is the medial axis and WIDTH is the width.
%   
%   SHAPE = TZ_MXS2MXP(MEDAXIS,WIDTH,PARAM) lets the user specify
%   parameters for splines. See TZ_MEDAXSPFEAT for more details.
%
%   See also TZ_MEDAXSPFEAT

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_mxs2mxp','ml_mxs2mxp'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

shape.format = 'mxp';

feat = tz_medaxspfeat(medaxis,width,param);

shape.length = feat{1};
shape.spmedaxis = tz_feat2sp(feat{2},feat{3});
shape.spwidth = tz_feat2sp(feat{4},feat{5});

