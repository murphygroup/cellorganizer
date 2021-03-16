function pts=ml_tracecontour(img,s)
%ML_TRACECONTOUR Trace contour of the object.
%   PTS = ML_TRACECONTOUR(IMG) returns a set of points that is the
%   contour of the object in image IMG. The priorities of tracing are
%   different for each direction:
%        4 3 2
%        5 0 1
%        6 7 8
%   The next direction also depends on current tracing direction.
%
%   PTS = ML_TRACECONTOUR(IMG,S) starts tracing the contour at point S.
%
%   See also

%   16-MAR-2005 Initial write  T. Zhao
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
    error('1 or 2 arguments are required')
end

if isempty(img)
    pts=[];
    return;
end

img = double(img);

if nargin<2
    s=[];
end

if isempty(s)
    [sx,sy]=find(img==1);
    s=[sx(1),sy(1)];
end

pts=s;
curpt=s;
nbstartpos = 1;

while 1
    img(curpt(1),curpt(2))=2;
%     imshow(img,[]);
%     drawnow
    nbsubs = [0,1;-1,1;-1,0;-1,-1;0,-1;1,-1;1,0;1,1];
    nbidx = ml_getnbidx(nbstartpos);
    nbsubs = nbsubs(nbidx,:);
    
    nbs=ml_addrow(nbsubs,curpt);
%     [curpt+[0,1];curpt+[-1,1];curpt+[-1,0];curpt+[-1,-1];...
%         curpt+[0,-1];curpt+[1,-1];curpt+[1,0];curpt+[1,1]];

    invalidIdx = find(nbs(:,1)<=0 | nbs(:,2)<=0 | ...
        nbs(:,1)>size(img,1) | nbs(:,2)>size(img,2));
    nbs(invalidIdx,:) = [];
    nbidx(invalidIdx) = [];
    
    for i=1:size(nbs,1)
        nbpt=img(nbs(i,1),nbs(i,2));
        found=0;
        if nbpt==1
            curpt=nbs(i,:);
            pts=[pts;curpt];
            nbstartpos = nbidx(i);
            found=1;
            break
        end
    end

    if found==0
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nbidx = ml_getnbidx(startpos)

reorder = [4 5 3 6 2 7 1 8];
idx = [6 7 8 1 2 3 4 5];
newidx = ml_shift(idx,1-startpos);
nbidx = newidx(reorder);
% 
% switch startpos
%     case 1
%         nbidx = [1 2 8 3 7 4 6 5];
%     case 2
%         nbidx = [2 3 1 4 8 5 7 6];
%     case 3
%         nbidx = [3 4 2 5 1 6 8 7];
%     case 4
%         nbidx = [4 5 3 6 2 7 1 8];
%     case 5
%         nbidx = [5 6 4 7 3 8 2 1];
% end
% 

