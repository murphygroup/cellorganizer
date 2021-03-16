function combfeats = ml_combccfeats(combcellcodes)
%ML_COMBCCFEATS Calculate features for cell shapes.
%   COMBFEATS = ML_COMBCCFEATS(COMBCELLCODES) returns the combined features
%   of cell codes COMBCELLCODES. Each cell code has following features:
%       'perimeter ratio'
%       'area ratio
%       'cell contour feats' 
%       'nuc eccentricity'
%       'cell-nuc angle' 
%       'cell-nuc distance along major axis' 
%       'cell-nuc distance along minor axis'
%       'cell major axis length' 
%       'cell minor axis length'
%       'nuc major axis length' 
%       'nuc minor axis length'
%   
%   See also

%   22-Apr-2005 Initial write  T. Zhao
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

 
if nargin < 1
    error('Exactly 1 argument is required')
end

for i=1:length(combcellcodes)
    s=combcellcodes{i};
    combfeats(i,1)=size(s.cellcontour,1)/size(s.nuccontour,1);
    combfeats(i,2)=s.cellarea/s.nucarea;
    combfeats(i,3:5)=ml_contourfeat(s.cellcontour);
    combfeats(i,6)=s.nucecc;
    combfeats(i,7)=mod(s.cellmangle-s.nucmangle,180);
    if (combfeats(i,7)>90)
        combfeats(i,7)=180-combfeats(i,7);
    end
    
    combfeats(i,7)=sign(s.cellmangle-s.nucmangle)*combfeats(i,7);
    
    diffcenter=s.nuccenter-s.cellcenter;
    
    theta=s.cellmangle*pi/180;
    mvec=[cos(theta),sin(theta)];
     
    combfeats(i,8)=abs(diffcenter*mvec');
    mproj=s.cellcontour*mvec';
    combfeats(i,10)=max(mproj)-min(mproj);
    
    mvec=[sin(theta),-cos(theta)];
    combfeats(i,9)=abs(diffcenter*mvec');
    mproj=s.cellcontour*mvec';
    combfeats(i,11)=max(mproj)-min(mproj);
    
    theta=s.nucmangle*pi/180;
    mvec=[cos(theta),sin(theta)];
    mproj=s.nuccontour*mvec';
    combfeats(i,12)=max(mproj)-min(mproj);
    
    mvec=[sin(theta),-cos(theta)];
    mproj=s.nuccontour*mvec';
    combfeats(i,13)=max(mproj)-min(mproj);
    
    combfeats(i,14)=combfeats(i,10)/combfeats(i,12);
end
