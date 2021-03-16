function h=tz_bwhand(img)
%TZ_BWHAND Handedness of an image. (Under modification)
%   H = TZ_BWHAND(IMG) returns the handedness of the binary version of 
%   an image IMG. There are three possible values: -1, 0, 1.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

img=img>0;

mom = tz_bwmoment(img);

center=[mom.cx,mom.cy];

theta = .5 * atan((mom.mu02 - mom.mu20)/2/mom.mu11)+sign(mom.mu11)*pi/4;%+pi/2;
[x,y]=find(img>0);

ntheta=[cos(theta),sin(theta)];
imgskew1=skewness([x,y]*ntheta');

ntheta=[cos(pi/2+theta),sin(pi/2+theta)];
imgskew2=skewness([x,y]*ntheta');

if any([imgskew1 imgskew2]==0)
    h=0;
    return
end

h=sign(imgskew1)*sign(imgskew2);