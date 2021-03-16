function [ratios, nucpts, cellpts] = ml_calccellnucratios3d(cellbody,nucbody,da)
%ML_CALCCELLNUCRATIOS3D 

% Created 28-Mar-2012 from ml_parsecell by T. Zhao
% Copyright (C) 2012 Murphy Lab
% Carnegie Mellon University
%
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

if nargin < 3
    error('CellOrganizer: Wrong number of input arguments. At least 3 arguments are required.');
end

len = sqrt((size(cellbody,1)^2+size(cellbody,2)^2));

%icaoberg march 29, 2012
%celledge=bwperim(cellbody,26);
%nucedge=bwperim(nucbody,26);
for i=1:1:size(cellbody,3)
    celledge(:,:,i)=bwperim(cellbody(:,:,i),8);
    nucedge(:,:,i)=bwperim(nucbody(:,:,i),8);
    celledge(:,:,i)= bwmorph( celledge(:,:,i),'dilate');
    nucedge(:,:,i)= bwmorph( nucedge(:,:,i),'dilate');
    celledge(:,:,i)= bwmorph( celledge(:,:,i),'thicken');
    nucedge(:,:,i)= bwmorph( nucedge(:,:,i),'thicken');
end

%murphy april 3, 2012
%this is done so that the cell body has top and botton non zero slices
celledge(:,:,1)=cellbody(:,:,1);
celledge(:,:,size(cellbody,3))=cellbody(:,:,size(cellbody,3));

%icaoberg march 29, 2012
nuccenter = regionprops( nucbody, 'Centroid' );
nuccenter = round( nuccenter.Centroid )

if ~isempty(celledge)
    discell=[]; disnuc=[]; nucpts=[]; cellpts=[];
    %for every zslice in the cell body
    for islice=1:size(cellbody,3)
      if (sum(sum(nucedge(:,:,islice)))>0)
        %islice
        %for every step given step size da
        for a=0:da:360-da
            %a
            %get points along a line from the nuclear center and a point given
            %by angle a and length len
            pts=ml_getlinept2(nuccenter(1:2),a,len);

%            a,size(pts)
%            if (islice==11 & a==226) keyboard, end

            %find points where that line intersects nuclear perimeter in the 
            %current slice and closest to nuclear center
            [point,distance]=ml_imgptspixelcell(nucedge(:,:,islice),pts,nuccenter(1:2));
            point=[point islice];
            distance=sqrt(sum((point-nuccenter).^2));
            disnuc=[disnuc;distance];
            nucpts=[nucpts;point];
            
%            % get points along a line between nuccenter and point
            pts=ml_getline2pts(nuccenter,point,len);
%            keyboard
%            point,nuccenter,size(pts)

            % find point along that line intersecting cell perimeter and
            % closest to nuclear center
            [point,distance]=ml_imgptspixel3(celledge,pts,nuccenter);
%            distance=sqrt(sum((point-nuccenter).^2));
            discell=[discell;distance];
            cellpts=[cellpts;point];
        end
      end
    end
end        

ratios=disnuc./discell;
%icaoberg april 2, 2012
ratios=ratios';
