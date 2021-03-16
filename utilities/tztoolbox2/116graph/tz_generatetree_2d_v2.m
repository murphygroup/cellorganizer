function [V,E]=tz_generatetree_2d_v2(nnode,dist_mean,dist_var,pdeg)
%TZ_GENERATETREE_2D_V2 Generate a minimal spanning tree.
%   V = TZ_GENERATETREE_2D_V2(NNODE,DIST_MEAN,DIST_VAR,PDEG) is similar
%   with TZ_GENERATETREE_2D. But it does not take angles as parameters.
%   It also uses fast neighbor searching algorithm.
%   
%   [V,E] = TZ_GENERATETREE_2D_V2(...) also returns edges.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('Exactly 4 arguments are required')
end

V(1,:)=[256,256];
distimg=zeros(512,512);
distimg(V(1),V(2))=1;
distimg=bwdist(distimg);

d=normrnd(dist_mean,dist_var);
[mindist,pos1,pos2]=tz_min(abs(distimg-d));
V(2,:)=[pos1,pos2];
mindist=distimg(pos1,pos2);
E(1,:)=[1,2,mindist];

degrees=[1 1];

if nnode<=2
    tz_plotgraph(V,E,[],'x','-',[])
    return;
end

degdistr=tz_mnornd(nnode-2,pdeg,1);

while any(degdistr(1:end-1)-degdistr(2:end)<0)
    degdistr=tz_mnornd(nnode-2,pdeg,1);
end

degs=[];
for i=1:length(degdistr)
    degs=[degs,ones(1,degdistr(i))+i];
end
degrsel=randperm(length(degs));

for i=3:nnode
    degsel=degrsel(i-2);
    inc=1;
    cadpt=find(degs(degsel)-degrees==1);
    
    while isempty(cadpt)
        degsel=degrsel(i-2+inc);
        inc=inc+1;
        cadpt=find(degs(degsel)-degrees==1);
    end
    
    degs(degsel)=0;
    
    succ=0;
    s=randperm(length(cadpt));
    k=1;
    ds=normrnd(dist_mean,dist_var);
    if ds<3
        ds=3;
    end
    
    %degrees
    while succ==0
        [i,ds]
        spt=V(cadpt(s(k)),:);
        distimg=zeros(512,512);
        distimg(spt(1),spt(2))=1;
        distimg=bwdist(distimg);
                
        [mdist,pos1,pos2]=tz_min(abs(distimg-ds));
        mindist=distimg(pos1,pos2);
        [cadsi,cadsj]=find(distimg==mindist);
        
        [index,ignore]=nn_search(V, nn_prepare(V), [cadsi,cadsj], 1);
        cadsi=cadsi(index==cadpt(s(k)));
        cadsj=cadsj(index==cadpt(s(k)));
        axis equal
        tz_plotgraph(V,E,[],'x','-',[])
        hold on
        plot(cadsi,cadsj,'.');
        
        drawnow
        hold off
        perm=randperm(length(cadsi));
         
        for m=1:length(cadsi)
            newv=[cadsi(perm(m)),cadsj(perm(m))];
            
            succ=tz_validatetredge(V,E,cadpt(s(k)),newv); 
            if succ==1
                k=k-1;
                break
            end
        end
        
        k=k+1;
        if k>length(s)
            k=1;
            ds=normrnd(dist_mean,dist_var);
            if ds<3
                ds=3;
            end
        end
    end
    V(i,:)=newv;
    E(i-1,:)=[cadpt(s(k)),i,mindist];
    %V
    degrees=[degrees,1];
    degrees(cadpt(s(k)))=degrees(cadpt(s(k)))+1;
    tz_plotgraph(V,E,[],'x','-',[])
    hold on
    drawnow
    
end



% figure
E=tz_mstree(V,'eu',1);
E=E';
% tz_plotgraph(V,E,[],'x','-',[])