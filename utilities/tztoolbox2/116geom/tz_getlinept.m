function pts=tz_getlinept(s,t)
%TZ_GETLINEPT Obsolete. See ML_GETLINEPT.
%   PTS = TZ_GETLINEPT(S,T) returns the coordinates of points on a line 
%   segment from S to T, which both are integer vectors with length 2.
%   PTS has two columns indicating X and Y coordinates.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_getlinept','ml_getlinept'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

if any(round([s t])~=[s t])
    error('the points must be integers');
end

if s(1)==t(1)
  if s(2)<t(2)
    pts=[zeros(t(2)-s(2)+1,1)+s(1),(s(2):t(2))'];
  end
  if s(2)==t(2)
    pts=s;
  end
  if s(2)>t(2)
    pts=[zeros(s(2)-t(2)+1,1)+s(1),(s(2):-1:t(2))'];
  end
  return;
end

if s(2)==t(2)
  pts=tz_getlinept([s(2),s(1)],[t(2),t(1)]);
  pts=[pts(:,2),pts(:,1)];
end

if all(s<t)
    
dx=t(1)-s(1);
dy=t(2)-s(2);

if dy < 0
  dy = -dy;
  stepy = -1;
else 
  stepy = 1;
end

if (dx < 0) 
  dx = -dx;  
  stepx = -1;
else 
  stepx = 1;
end

dy=dy*2;
dx=dx*2;

pts(1,:)=s;


if (dx > dy)
  fraction = dy - dx/2;
  while (s(1) ~= t(1))
    if fraction >= 0
      s(2)=s(2)+stepy;
      fraction=fraction-dx;
    end
    s(1)=s(1)+stepx;
    fraction=fraction+dy;
    pts=[pts;s];
  end
else
  fraction= dx - dy/2;
  while (s(2) ~= t(2))
    if (fraction >= 0)
      s(1)=s(1)+stepx;
      fraction=fraction-dy;
    end
    s(2)=s(2)+stepy;
    fraction=fraction+dx;
    pts=[pts;s];
  end
end

end

if s(1)>t(1) & s(2)<t(2)
    pts=tz_getlinept([t(1)-s(1),s(2)],[0,t(2)]);
    pts(:,1)=t(1)-pts(:,1);
end

if all(s>t) |  (s(1)<t(1) & s(2)>t(2))
  pts=tz_getlinept(t,s);
  pts=flipud(pts);
end





