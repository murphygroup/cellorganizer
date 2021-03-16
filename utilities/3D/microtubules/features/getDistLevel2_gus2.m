function [l_img] = getDistLevel2_gus(image,segcell, segdna ,numlevels,check)

% Copyright (C) 2015-2016 Murphy Lab
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

%What is this function supposed to do? grj 10/4/14
if ~exist('check','var')
    check = 1;
end

if isa( image, 'double' )
    segcell = double(segcell);
end

if isa( image, 'uint8' )
    segcell = uint8(segcell);
end
image = segcell.*image;

[D,L] = bwdist(~segcell);
%unique_D = unique(D);
%unique_D(1) = [];

l_img = segcell*0;

ww = find(segdna>0);
D(ww) = -1;
m_level = linspace(1,max(D(:)),numlevels+1);

for i=1:numlevels
    ww = find( D >= m_level(i) & D<m_level(i+1) );
    l_img(ww) = i;
end

