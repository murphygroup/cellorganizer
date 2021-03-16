function f = tz_rdistpca(combcellcodes,param)
%TZ_RDISTPCA Train the the PCA model for distance ratios.
%   F = TZ_RDISTPCA(COMBCELLCODES) returns the trained PCA model for the
%   input cell array of cell codes COMBCELLCODES. F is a structure from
%   TZ_TRAINLK.
%   
%   F = TZ_RDISTPCA(COMBCELLCODES,PARAM) specifies parameters for training.
%   PARAM is a structure that has the following fields:
%       'startangle' - a string that determines how to align the ratio
%           vector. 'cell' means the major angle of the cell and 'nuc'
%           means the major angle of the nucleus.
%       'tz_trainlk' - parameters for the function TZ_TRAINLK
%
%   See also

%   10-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

warning(tz_genmsg('of','tz_rdistpca','ml_rdistpca'));

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('startangle','cell','tz_trainlk',struct([])));

for i=1:length(combcellcodes)
    cellcode = combcellcodes{i};
    rdist = cellcode.celldist./cellcode.nucdist;
    switch param.startangle
        case 'cell'
            startAngle = cellcode.cellmangle;
        case 'nuc'
            startAngle = cellcode.nucmangle;
    end
    normrdist(i,:) = tz_shiftdist(rdist,startAngle);
end

% f = tz_trainlk(normrdist,param.tz_trainlk);
f = ml_estpdf(normdist,param.ml_estpdf);
