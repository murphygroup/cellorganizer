function image_features=tz_calfeats_3d(images,feat_fullpath)
%TZ_CALFEATS_3D Obsolete.

%function tz_calfeats_3d(images,feat_fullpath)
%OVERVIEW:
%   Calculate 3d features of the image stored in images.mat, then save the result in feat_fullpath
%PARAMETERS:
%   images - matlab file storing 3d matrix
%   feat_fullpath - feature saved path
%RETURN:
%   image_features - calculated features
%DESCRIPTION
%   Calculate 14 morphological features (3D-SLF 9.1~9.8, 9.15~9.20)
%
%HISTORY
%   ??-???-???? Initial write TINGZ
%   14-DEC-2003 Modified TINGZ
%   19-MAY-2004 Modified TINGZ
%       - add return value

error(tz_genmsg('of','tz_calfeats_3d'));

% Compute features
image_features = ml_3dfeatures( images, [], [], [0.1075 0.1075 0.5]);

% Save features for this image
save( feat_fullpath,'image_features');
