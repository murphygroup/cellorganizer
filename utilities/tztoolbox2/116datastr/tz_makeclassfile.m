function tz_makeclassfile(filename,features)
%TZ_MAKECLASSFILE Generate a file for classification.
%   TZ_MAKECLASSFILE(FILENAME,FEATURES) generates a matlab file FILENAME
%   for the cell array of feature matrices FEATURES. This is previously 
%   designed for HK_EXPRBFSVM. Now it is almost useless.

%   ??-???-???? Initial write T. Zhao
%   26-MAR-2003 Modified T. Zhao
%   14-DEC-2003 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

classnum=length(features);
all_features={};

for i=1:classnum
    all_features=tz_addclass(all_features,features{i},-1);
end

save(filename,'all_features');