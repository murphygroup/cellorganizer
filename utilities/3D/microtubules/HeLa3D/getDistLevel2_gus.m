function [l_img] = getDistLevel2_gus(image,imnum,numlevels,check)

if ~exist('check','var')
   check = 1;
end

[protim3,Dbothfin,segdna,segcell,dnaim3,cellim3,imgcent_coordinate] = getrealimage_hela(imnum);
if check && ismember(imnum,[1,6,14,19,31,52])  %%
   segcell = cat(3,segcell,zeros(size(segcell,1),size(segcell,2),16-size(segcell,3)));
end
image = segcell.*image;

[D,L] = bwdist(~segcell);
%unique_D = unique(D);
%unique_D(1) = [];

l_img = segcell*0;

ww = find(segdna>0);
D(ww) = -1;
m_level = linspace(1,max(D(:)),numlevels+1);

for i=1:numlevels
   
    ww = find( D >= m_level(i) & D<m_level(i+1) );
    l_img(ww) = i;
end


%interval = floor(length(unique_D)/numlevels);

%levels = cell(numlevels,1);

%for i = 1:numlevels-1
%    [idxr,idxc,idxz] = ind2sub(size(D),find(D>unique_D(1+interval*(i-1))-1e-6 & D<unique_D(interval*i+1)-1e-6));
%    idx = find(D>unique_D(1+interval*(i-1))-1e-6 & D<unique_D(interval*i+1)-1e-6);
%    l_img(idx) = i;
    
%    levels{i} = [idxr(:),idxc(:),idxz(:),image(idx)];
%end
%    i = numlevels;
%    [idxr,idxc,idxz] = ind2sub(size(D),find(D>unique_D(1+interval*(i-1))-1e-6 & D<unique_D(end)+1e-6));
%    idx = find(D>unique_D(1+interval*(i-1))-1e-6 & D<unique_D(end)+1e-6);
%    levels{i} = [idxr(:),idxc(:),idxz(:),image(idx)];
%    l_img(idx) = i;