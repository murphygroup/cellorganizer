function features = ml_obj2feat( protobj, DNACOF, ...
				 protimg,  ...
				 protbin, dnabin, ...
				 scale, scaledist)

% FEATURES = ML_OBJ2FEAT( PROTOBJ, DNACOF, ...
%                         PROTIMG, ...
%                         PROTBIN, DNABIN, SCALE)
%
% Calculates features for a cell based on a list of objects in the
% cell. The list of object should be one like returned by
% ml_3dfindobj, which takes a 3D image and finds objects in it.
% PROTOBJ is the list of objects in the protein image.
% PROTIMG is the clean image from which the
% objects were found. This is necessary in order to calculate
% the COF accurately, because PROTOBJ does not include pixel
% values.
% PROTBIN, DNABIN and CELLBIN are the respective binary images that 
% were used to find the objects in the images.
%
% The features calculated are (in the order they are returned):
%
% 1) No. of objects in the image
% 2) Euler number of the image
% 3) Average object volume (average number of above threshold voxels per object)
% 4) Std.Dev of object volumes
% 5) Ratio of Max(object volumes) to Min(object volumes)
% 6) Average Obj to protein COF distance
% 7) Std.Dev of Obj to protein COF distances
% 8) Ratio of Max to Min obj distance from prot COF
% 9) Average Obj to DNA COF distance
%10) Std.Dev of Obj to DNA COF distances
%11) Ratio of Max to Min obj distance from DNA COF
%12) Distance between Prot COF and DNA COF
%13) Ratio of protein volume to DNA volume
%14) Fraction of protein fluorescence overlapping with DNA fluorescence
%15-28) horizontal-vertical directional features based on 6-12 above

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

% Modified from ml_obj2feat to check whether DNA chanel is available.
% Xiang Chen Nov 18, 2004

if( ~exist( 'scaledist', 'var'))
     if( ndims(protimg)== 3)
         scaledist = [1 1 203/48.8];
     else
         scaledist = ones(1,ndims(protimg));
     end
end   % Feb 07, YH, added the default scale (3D Hela scaledist)

NObj = length( protobj);
features(1) = NObj; % Number of Objects in the image
holes = [];
sizes = [];
for o = 1 : NObj
    holes = [holes protobj{o}.n_holes];
    sizes = [sizes protobj{o}.size];
end
total_holes = sum( holes);
features(2) = NObj - total_holes; % Euler Number of the image
features(3) = mean( sizes)*prod(scale); % Average object volume
features(4) = std( sizes)*prod(scale); % Standard Deviation of obj volumes

if( min( sizes) > 0)
    features(5) = max( sizes) / min( sizes); % Ratio of Max volume to Min volume
else
    error('Why is min(objsize) = 0 ?');
end

% find COF stuff
[ProtCOF, ProtCOFs] = ml_findCOFs( protobj, protimg);
ProtCOFDists = ml_eucdist( ProtCOF, ProtCOFs, scaledist); % this and subsequent 'scaledist' added, jnewberg 11/24/05

% COF features relating object distances to Protein COF
[features(6),features(7),features(8)] = ml_ObjCOFfeats( ProtCOFDists);


% COF features relating Protein to DNA
if (~isempty(DNACOF))
    % Average dist, stddev dist and max/min dist.
    ProtDNADists = ml_eucdist( DNACOF, ProtCOFs, scaledist);
    [features(9),features(10),features(11)] = ml_ObjCOFfeats( ProtDNADists);
    % Distance between Prot COF and DNA COF
    features(12) = ml_eucdist( ProtCOF, DNACOF, scaledist);
    % Ratio of protein volume to DNA volume
    % Find protein volume
    protvolume = length(find(protbin));
    dnavolume = length(find(dnabin));
    features(13) = protvolume / dnavolume;
    % Fraction of protein fluorescence overlapping with DNA fluorescence
    features(14) = length(find(protbin & dnabin)) / protvolume;
else
    features(9:14) = nan;
end

%%%%%%%%%%%% 3D vertical-horizontal sensitive features%%%%%%%%%%%%%%%%%%
% Get the horizontal and vertical components of distances of
% objects from COF of protein and DNA
ProtCOFDistsH = ml_eucdist( ProtCOF(1:2,:), ProtCOFs(1:2,:),scale(:,1:2));
ProtCOFDistsV = ml_eucdist( ProtCOF(3,:), ProtCOFs(3,:),scale(3));
if (~isempty(DNACOF))
    ProtDNADistsH = ml_eucdist( DNACOF(1:2,:), ProtCOFs(1:2,:),scale(:,1:2));
    ProtDNADistsV = ml_eucdist( DNACOF(3,:), ProtCOFs(3,:),scale(3));
end
% COF features relating object distances to Protein COF
[features(15),features(16),features(17)] = ml_ObjCOFfeats( ProtCOFDistsH);
[features(18),features(19),features(20)] = ml_ObjCOFfeats( ProtCOFDistsV);

if (~isempty(DNACOF))
    % COF features relating Protein to DNA
    % Average dist, stddev dist and max/min dist.
    [features(21),features(22),features(23)] = ml_ObjCOFfeats( ProtDNADistsH);
    [features(24),features(25),features(26)] = ml_ObjCOFfeats( ProtDNADistsV);
    % Distance between Prot COF and DNA COF
    features(27) = ml_eucdist( ProtCOF(1:2,:), DNACOF(1:2,:),scale(:,1:2)); % horizontal
    features(28) = (ProtCOF(3,:) - DNACOF(3,:))*scale(3); % vertical
else
    features(21:28) = nan;
end
