function [V,E]=tz_generatetree_2d_v3(nnode,dist_mean,dist_var,ndeg,nmbrch)
%TZ_GENERATETREE_2D_V2 Generate a minimal spanning tree.
%   V = TZ_GENERATETREE_2D_V2(NNODE,DIST_MEAN,DIST_VAR,NDEG,NMBRCH) is 
%   similar with TZ_GENERATETREE_2D. But it does not take angles as 
%   parameters. It also uses fast neighbor searching algorithm. NMBRACH
%   is the number of nodes in the main branch. NDEG(I) is the number of 
%   degrees no less than I+1.
%   
%   [V,E] = TZ_GENERATETREE_2D_V2(...) also returns edges.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 5
    error('Exactly 5 arguments are required')
end  

%The first node
V(1,:)=[256,256];

%create distance field
distimg=zeros(512,512);
distimg(V(1),V(2))=1;
distimg=bwdist(distimg);

%generate the first edge and the second node
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

degdistr=ndeg;

degs=[];
for i=1:length(degdistr)
    degs=[degs,ones(1,degdistr(i))+i];
end



maxiter=1000;
iter=maxiter+1;

while iter>maxiter
    iter=1;
    if nmbrch>2
        degrsel=[1:nmbrch-2,randperm(sum(degdistr)-nmbrch+2)+nmbrch-2];
    else
        degrsel=randperm(length(degs));
    end
    for i=3:nnode
        rsel=find(degrsel>0);
        degsel=degrsel(rsel(1));
        inc=1;
        cadpt=find(degs(degsel)-degrees==1);
        
        while isempty(cadpt)
            
            degsel=degrsel(rsel(inc+1));
            
            inc=inc+1;
            cadpt=find(degs(degsel)-degrees==1);
        end
        
        degs(degsel)=0;
        degrsel(rsel(inc))=0;
        
        succ=0;
        s=randperm(length(cadpt));
        k=1;
        ds=normrnd(dist_mean,dist_var);
        if ds<3
            ds=3;
        end
        
        %degrees
        while succ==0
            iter=iter+1;
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
            axis equal
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
            if iter>maxiter
                break;
            end
            
        end
        if iter>maxiter
            break;
        end
        V(i,:)=newv;
        E(i-1,:)=[cadpt(s(k)),i,mindist];
        %V
        degrees=[degrees,1];
        degrees(cadpt(s(k)))=degrees(cadpt(s(k)))+1;
        tz_plotgraph(V,E,[],'x','-',[])
        hold on
        drawnow
        
        if nmbrch>2
            if size(V,1)>nmbrch
                degrees(degrees==1)=0;
                nmbrch=0;
            end
        end
    end
    
end

figure
E=tz_mstree(V,'eu',1);
E=E';
tz_plotgraph(V,E,[],'x','-',[])