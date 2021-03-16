function dist=tz_wnddist(wnd1,wnd2,mask,option)

%function dist=tz_wnddist(wnd1,wnd2)

if sum(mask(:))==0
    dist=Inf;
    return;
end

wndsize=prod(size(wnd1));

x1=reshape(wnd1,[1,wndsize]);
x2=reshape(wnd2,[1,wndsize]);
mask=reshape(mask,[1,wndsize]);
x1(mask==0)=[];
x2(mask==0)=[];

switch(option)
case 'pixel'
    dist=sum(sum((x1-x2).^2));
case 'hist'   
    x=[x1,-x2];
    dist=sum((kdest1d(x1,x','knorm',0)-kdest1d(x2,x','knorm',0)).^2);
case 'nssd'
    dist=mean((x1-x2).^2);
end