function nucsurf = nucimg2nucsurf( nucimg )
%NUCIMG2NUCSURF Helper method that converts a nucleus image to a nucleus
%surface representation
%
%Input: nucleus image (3D matrix)
%Output: nucleus surface representation (2D matrix where each row is the polar coordinate index of the surface)

% Author: Yue Yu(yuey1@andrew.cmu.edu)
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% October 12, 2012 I. Cao-Berg Documented method, added check of input
% arguments, encapsulated method in try/catch method so that method returns
% empty array if it fails
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

%icaoberg 10/12/2012
%Yue Yu 10/12/2012 fixed a interpolation bug 
if nargin ~= 1
    warning( 'CellOrganizer: Wrong number of input arguments' );
    nucsurf = [];
elseif isempty( nucimg )
    warning( 'CellOrganizer: Input argument cannot be empty' );
    nucsurf = [];
elseif size( nucimg, 3 ) <= 1
    warning( 'CellOrganizer: Input argument must be a 3D array' );
    nucsurf = [];
else
    %icaoberg 10/12/2012
    try
        %yuey1 10/10/2012
        [r,c,z] = size(nucimg);
        
        nucsurf = zeros(z,361);
        for i = 1 : z
            tmp_slice = bwperim(nucimg(:,:,i));
            tmp_ind = find(tmp_slice ~= 0);
            if isempty(tmp_ind)
                continue;
            end
            [yy,xx] = ind2sub([r,c],tmp_ind);
            yy = yy - c/2;
            xx = xx - c/2;
            [TH,R] = cart2pol(yy,xx);
            [~,tmpi] = sort(TH);
            R = R(tmpi);
            if length(R) >= 361
                rr = interp1(1:length(R),R,1:length(R)/361:length(R));
            else
                rr = interp1(1:length(R),R,1:length(R)/362:length(R));
            end
            nucsurf(i,:) = fliplr(circshift(rr',270)');
        end
    catch
       warning( 'CellOrganizer: Unable to make nuclear surface representation for the given image' );
       nucsurf = [];
    end
end