function speed = tz_constflow(img1,img2,wnd)
%TZ_CONSTFLOW 
%   SPEED = TZ_CONSTFLOW(IMG1,IMG2)
%   
%   SPEED = TZ_CONSTFLOW(IMG1,IMG2,WND)
%   
%   See also

%   31-Oct-2005 Initial write T. Zhao

if nargin < 2
    error('2 or 3 arguments are required')
end

if nargin < 3
    wnd = ones(size(img1));
end

img1 = double(img1);
img2 = double(img2);

[Iy,Ix] = gradient(img1);
Ix = Ix(find(wnd==1));
Iy = Iy(find(wnd==1));
A = [Ix,Iy];
speed = [0;0];
MINERROR = 1e-3;

It = tz_imsub(img2,img1);
It = It(find(wnd==1));
speed2 = -inv(A'*A)*A'*It;

for k=1:30
    if all(abs(speed2)<MINERROR)
        break;
    end
    speed = speed+speed2
    
    img2 = tz_imtranslate(img2,-speed2','bilinear');
    It = tz_imsub(img2,img1);
    It = It(find(wnd==1));
    speed2 = -inv(A'*A)*A'*It;
    
end