function [pt,distnce]=ml_imgptspixel3(img,pts,center)
%ML_IMGPTSPIXEL3 Get pixel values at specified points.
%   PS = ML_IMGPTSPIXEL3(IMG,PTS) returns a vector of pixel values from
%   the image IMG. PS(I) is the gray level of the pixel at position
%   [PTS(I,1),PTS(I,2)] in IMG.

% March 28, 2012 R.F.Murphy created from ml_imgptspixel by T. Zhao
% April 2, 2012 I. Cao-Berg Modified the method so that if it doesn't get a hit
%                           from the coordinates in pts, it looks for the closest
%                           coordinate on the img to the end of the line defined by pts 
%
% Copyright (C) 2012  Murphy Lab
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

if nargin < 2
    error('Exactly 2 arguments are required')
end

extra=[zeros(size(pts,1),1) zeros(size(pts,1),1) ones(size(pts,1),1)];
ptscopy=[pts-extra; pts; pts+extra];

imgsize=size(img);
ptscopy(ptscopy(:,1)<=0 | ptscopy(:,1)>imgsize(1),:)=[];
ptscopy(ptscopy(:,2)<=0 | ptscopy(:,2)>imgsize(2),:)=[];
ptscopy(ptscopy(:,3)<=0 | ptscopy(:,3)>imgsize(3),:)=[];

if isempty(ptscopy)
    error('empty points');
end

idx=sub2ind(imgsize,ptscopy(:,1),ptscopy(:,2),ptscopy(:,3));

ps=img(idx);
int3=find(ps>0);                      

%keyboard

if ~isempty( int3 )
  dis3=sqrt(sum((ptscopy(int3,:)-repmat(center,size(int3))).^2,2));
  distnce=min(dis3);
  pt = ptscopy( int3(find( dis3 == distnce)),:,: )';
  pt = pt(:,1)'; % in case there is more than one point that is the same dist

%  imshow(img(:,:,pt(3))); hold on;
%  plot(ptscopy(:,1),ptscopy(:,2),'r-');
%  plot(ptscopy(int3,1),ptscopy(int3,2),'g+'); hold off; pause(0.1);


else
    disp('ml_imgptspixel3: No intersection.');
for i=1:size(img,3)
    imshow(img(:,:,i)); hold on;
    xx=find(ptscopy(:,3)==i);
    if (length(xx)>0) plot(ptscopy(xx,1),ptscopy(xx,2),'r-'); end
    pause;
    hold off;
end
   %since there are not hits, find the closest point from the end in ptscopy to the img
   %pt = findClosest( img, ptscopy );
end
