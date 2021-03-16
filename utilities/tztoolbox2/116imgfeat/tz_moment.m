function Out=tz_moment(img)
%TZ_MOMENT Central moments of an image.
%   TZ_MOMENT(IMG) returns the central moments of an image.
%   The return value is a structure containing the first, second
%   moments and the zero,first,second,third central moments.

%   ??-???-???? Modified from aw_moment TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add comments

total_intensity = sum(img(:));

if total_intensity>0
    intensity=img/total_intensity;
    posy=repmat(1:size(img,2),size(img,1),1);
    posx=repmat(1:size(img,1),size(img,2),1)';
    
    fxyz = sum(intensity(:));
    ix=posx.*intensity;
    iy=posy.*intensity;
    xfxyz = sum(ix(:));
    yfxyz = sum(iy(:));
    ix2=posx.*ix;
    iy2=posy.*iy;
    x2fxyz = sum(ix2(:));
    y2fxyz = sum(iy2(:));
    x3fxyz = sum(sum(posx.*ix2));
    y3fxyz = sum(sum(posy.*iy2));
    xyfxyz = sum(sum(posx.*iy));
    x2yfxyz = sum(sum(posy.*ix2));
    xy2fxyz = sum(sum(posx.*iy2));

    % Central Moments!
    
    xbar = xfxyz/fxyz;
    ybar = yfxyz/fxyz;
else
    xbar = 0;
    ybar = 0;
end

mu00 = fxyz;

mu20 = x2fxyz - (mu00 * xbar^2);
  
mu02 = y2fxyz - (mu00 * ybar^2);

mu11 = xyfxyz - (mu00 * xbar * ybar);

mu30 = x3fxyz - (3 * x2fxyz * xbar) + (2 * mu00 * xbar^3);
mu03 = y3fxyz - (3 * y2fxyz * ybar) + (2 * mu00 * ybar^3);

mu21 = x2yfxyz - (x2fxyz * ybar) - (2 * xyfxyz * xbar) + (2 * mu00 * xbar^2 * ybar);

mu12 = xy2fxyz - (y2fxyz * xbar) - (2 * xyfxyz * ybar) + (2 * mu00 * xbar * ybar^2);

%Output!


Out = struct('x2',x2fxyz,'y2',y2fxyz,'xbar',xbar,'ybar',ybar,'mu00',mu00,'mu11',mu11,'mu20', mu20,'mu02',mu02,'mu30',mu30,'mu03',mu03,'mu12',mu12,'mu21',mu21);

