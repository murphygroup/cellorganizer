function scale = tz_getresolution(dataset)
%TZ_GETRESOLUTION Resolution of microscopy image data set.
%   SCALE = TZ_GETRESOLUTION(DATASET) returns the resolution of image data
%   set specified by DATASET:
%       '3dhela' - 3D HeLa dataset
%       '3d3T3' - 3D 3T3 dataset
%   The unit of SCALE is micrometer.

%   20-Jul-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

switch(dataset)
    case '2dhela'
        scale = [0.23 0.23];
    case '3dhela'
        scale=[0.049 0.049 0.203];
    case '3d3T3'
        scale=[0.1075 0.1075 0.5];
    otherwise
        error('unknown dataset');
end
