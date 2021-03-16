function cpixels=tz_getknnpixels(img,mask,wnd,qmask,df,maxmasked,ep)

%function cpixels=tz_getknnpixels(img,mask,wnd,df,ep)

if isempty(mask)
    mask=ones(size(img));
end

k=floor(size(wnd,1)/2);

extimg=tz_extendimg(img,k,0);
extmask=tz_extendimg(mask,k,0);
cpixels=[];

for i=1:size(img,1)
    for j=1:size(img,2)
        maskwnd=tz_getnbwnd(extmask,[i,j],k,0);
        if ~isempty(qmask)
            qmask(k+1,k+1)=1;
            maskwnd=maskwnd.*qmask;
        end
        
        if maskwnd(k+1,k+1)==1
            if sum(sum(maskwnd==0))<=maxmasked
                imgwnd=tz_getnbwnd(extimg,[i,j],k,0);
                maskwnd(k+1,k+1)=0;
                if tz_wnddist(wnd,imgwnd,maskwnd,df)<=ep
                    cpixels=[cpixels imgwnd(k+1,k+1)];
                end
            end
        end
    end
end

