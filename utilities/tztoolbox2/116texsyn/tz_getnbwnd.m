function wnd=tz_getnbwnd(img,pos,k,prep)

%function wnd=tz_getnbwnd(img,pos,k)

if prep==1
    img=tz_extendimg(img,k);
end

wndpos=[pos(1):pos(1)+2*k;pos(2):pos(2)+2*k];
   
wnd=img(wndpos(1,:),wndpos(2,:)); 

