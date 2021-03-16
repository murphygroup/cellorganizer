function theta=ml_majorangle(img)

%ML_MAJORANGLE calculates the major angle an image.
%   THETA=ML_MAJORANGLE(IMG) returns the major angle of 
%   IMG with the unit radius. The angle is calculated
%   as the couter clockwise angle between x axis and 
%   the major axis. In the coodinate system using here,
%   x and y are row and column repectively.

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

%   01-NOV-2004 Modified T. Zhao
%   24-Mar-2005 Modified T. Zhao
%       - debugged
%   21-Feb-2006 Modified T. Zhao
%       - consider zero mu11

if min(img(:))==max(img(:))
    theta=0;
    warning('constant image');
end

%Calculate necessary moments
mom = ml_moment2(img);
weights=img(find(img>0));

center=[mom.cx,mom.cy];

if mom.mu11==0
    theta = 0;
else
    theta = .5 * atan((mom.mu02 - mom.mu20)/2/mom.mu11)+sign(mom.mu11)*pi/4;

    ntheta=[cos(theta),sin(theta)];
    [x,y]=find(img>0);
    x=x-mom.cx;
    y=y-mom.cy;

    %do projection and estimate skewness
    imgskew=ml_wmoment([x,y]*ntheta',weights,3);

    if imgskew<0
        theta=theta+pi;
    end
end
