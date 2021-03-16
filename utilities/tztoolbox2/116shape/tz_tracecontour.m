function pts=tz_tracecontour(img,s)
%TZ_TRACECONTOUR Obsolete. See ML_TRACECONTOUR.
%   PTS = TZ_TRACECONTOUR(IMG) returns a set of points that is the
%   contour of the object in image IMG. The priorities of tracing are
%   different for each direction:
%        4 3 2
%        5 0 1
%        6 7 8
%   The next direction also depends on current tracing direction.
%
%   PTS = TZ_TRACECONTOUR(IMG,S) starts tracing the contour at point S.
%
%   See also

%   16-MAR-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_tracecountour','ml_tracecontour'));

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
    nbidx = tz_getnbidx(nbstartpos);
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
function nbidx = tz_getnbidx(startpos)

reorder = [4 5 3 6 2 7 1 8];
idx = [6 7 8 1 2 3 4 5];
newidx = shift(idx,1-startpos);
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

