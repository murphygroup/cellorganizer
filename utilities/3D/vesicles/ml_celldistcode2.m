function [distcodes,coords,angles] = ml_celldistcode2( param, loc )
%ML_CELLDISTCODE2 Coding pixels by their distaces to cell and nuclues edge.
%   DISTCODES = ML_CELLDISTCODE2(param.cell,param.nucleus,OBJECTS,ISBIN) returs a
%   3-column matrix DISTCODES. Each row of DISTCODES has the form
%   [[nucleus distance,cell distance,gray level]. param.cell and param.nucleus are
%   [image]s for cell edge and nuclues edge. 
%   
%   ML_CELLDISTCODE2(param.cell,param.nucleus,OBJECTS,ISBIN,LOC) specifies what
%   regions are coded:
%       'cyt' - cytoplasm
%       'nuc' - nucleus
%       'all' - whole cell
%
%   [DISTCODES,COORDS,ANGLES] = ML_CELLDISTCODE(...) also returns the 
%   coordinates.

% March 19, 2012 R.F. Murphy and Devin Sullivan 
%   Created from  ml_celldistcode.m by T. Zhao
% grj 9/25/13 - changed from find(codemask==1)
%
% Copyright (C) 2007-2012 Murphy Lab
% Carnegie Mellon University
%
% August 4, 2012 D. Sullivan Passed param.nucleus and param.cell from parent call
%                            to save memory
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

if ~exist('loc','var')
    loc = 'cyt';
end

if any(size(param.cell)~=size(param.nucleus))
    error('Sizes of edge images do not match.');
end

switch loc
    case 'cyt'
        codemask = double(param.cell>0)-double(param.nucleus>0);
    case 'nuc'
        %devins ??/??/2012
        codemask = double(param.nucleus - bwperim(param.nucleus));
    case 'all'
        codemask = double(param.cell- bwperim(param.cell));
    otherwise
        error(['unrecognized location name: ' loc]);
end

if param.verbose
    disp( 'Calculating distances' );
end
ind = find(codemask>0); %grj 9/25/13 - changed from find(codemask==1)
[X,Y,Z] = ind2sub(size(codemask),ind);
coords = [X,Y,Z];

distcodes = [param.nucleardist(ind),param.celldist(ind),param.cell(ind)];

%Calculate angles
if param.verbose
   disp( 'Calculating angles' );
end

if nargout==3
    %D. Sullivan 7/17/13
    [center,mangle] = ml_objcenter(param.nucleus>0);
%     [center,mangle] = ml_objcenter(param.nucleus);
    %pt1 = center+[cos(mangle),sin(mangle)];
    %pt2 = center;
    [th,phi,r] = cart2sph(coords(:,1)-center(1),coords(:,2)-center(2),...
        coords(:,3)-center(3));
    th = th - mangle;
    angles.theta = th;
    angles.phi = phi;
    %angles = ml_normangle(ml_normangle( ...
    %   ml_normangle(th,'r2a'),'360'),'a2r');
    
    %icaoberg ??/??/2010
    %Old algorithm (slow)
    %for i=1:size(coords,1)
    %    angles2(i) = tz_findangle_2d(pt1,pt2,coords(i,:));
    %end
    %sum(abs(angles-angles2'))
end
