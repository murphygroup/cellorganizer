function SmallCellData=Process_SmallCell(img,loc,width,height)
position=[loc(1)-width/2,loc(2)-height/2,width,height];
%border processing
if position(1)<1
    position(1)=1;
end
if position(1)+width>size(img,2)
    position(1)=size(img,2)-width;
end
if position(2)<1
    position(2)=1;
end
if position(2)+height>size(img,1)
    position(2)=size(img,1)-height;
end
tmp=imcrop(img,position);
BW=im2bw(tmp,graythresh(tmp));
% Get the maximum concatenation field other set to 0
L = bwlabel(BW);
stats = regionprops(L);
Ar = cat(1, stats.Area);
ind = find(Ar ==max(Ar));
if(length(ind)==1)
    BW(L~=ind)=0;
    MAX_AR=max(Ar);
else
    MAX_AR=0;
end

% If the maximum connected area in the area to be pulled is greater than 0.3,note that this target is valid
if MAX_AR/width*height>0.3
    BW=medfilt2(BW,[2 2]);
    BW = imfill(BW, 'holes');
    Edge=edge(BW,'canny');
    [m,n]=find(Edge==1);
    if max(m) - min(m) < 3 || max(n) - min(n) < 3
        SmallCellData=struct();
        return;
    end
    m=m+loc(2)-height/2;n=n+loc(1)-width/2;
    SmallCellData=Ellipse_Fitting2([n,m]);
else
    SmallCellData=struct();
end

    