function L = seededwatershed(A,seeds,conn)
% L = SEEDEDWATERSHED(A,SEEDS,CONN)
% L is the output, watershed image regions (type double)
% A is the input gray level image (type uint16)
% SEEDS is the seed channel, a binary image)
% CONN is the connectivity, by default 4 for 2D, 6 for 3D
% 

% Copyright (C) 2008  Murphy Lab
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

if ~exist( 'conn','var')
    if ndims(A)==2
        conn = 4;
    elseif ndims(A)==3
        conn = 6;
    else
        error( 'Improper image dimensions');
    end
end

if size(A)~=size(seeds)
    error( 'Inputs are not compatible');
end

M = bwlabeln( seeds,conn);
L = watershed_meyer(A,conn,M);

% This is faster but does not work with touching seeds
% idx = unique(L(seeds==max(seeds(:))));
% L = ismember(L,idx);

% Slower than the above method, faster than below, and gives same as below
idx = find(seeds>0 & seeds<max(seeds(:)));
idx = unique(L(idx));
M = ismember(L,idx);
L=L>0;
L(M==1)=0;

% [r c] = find(seeds>0 & seeds<max(seeds(:)));
% L = ~imfill(~L,[r c],conn);
