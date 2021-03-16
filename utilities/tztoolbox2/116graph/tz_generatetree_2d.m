function [V,E]=tz_generatetree_2d(nnode,angle_mean,angle_var,...
    dist_mean,dist_var,pdeg)
%TZ_GENERATETREE_2D Generate a minimal spanning tree.
%   V = 
%   TZ_GENERATETREE_2D(NNODE,ANGLE_MEAN,ANGLE_VAR,DIST_MEAN,DIST_VAR,PDEG)
%   returns a two-column matrix of node coordinates for the generated tree.
%   The generated tree will have NNODE nodes. ANGLE_MEAN and ANGLE_VAR
%   determine the parameters of generating angles of the branches by mean
%   value and variance. DIST_MEAN and DIST_VAR determine the parameters
%   of generating lengths of the edges by mean and variance. PDEG is a
%   vector specifying the frequencies of node degrees. PDEG(I) is the 
%   estimated probability of degrees no less than I+1.
%
%   [V,E] = TZ_GENERATETREE_2D(...) also return edges of the tree.


%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 6
    error('Exactly 6 arguments are required')
end

V(1,:)=[0,0];
dist=normrnd(dist_mean,dist_var);
V(2,:)=[dist,0];
E(1,:)=[1 2 dist];
degrees=1;

orgV=V;
orgE=E;
maxiter=100000;

while(size(V,1)<nnode)
    succ=0;
    maxdeg=max(degrees);
    deg=find(tz_mnornd(1,pdeg,1)==1)+1;
    while all(deg-degrees~=1)
        deg=find(tz_mnornd(1,pdeg,1)==1)+1;
    end
    iter=1;
    
    while succ==0
        if iter>maxiter
            break;
        end
        cadpt=find(deg-degrees==1);
        s=randperm(size(cadpt,1));
        s=cadpt(s(1))+1;
        if isempty(angle_mean)
            a=unifrnd(0,2*pi);
        else        
            a=normrnd(angle_mean,angle_var);
        end
        dist=normrnd(dist_mean,dist_var);
        
        pts=tz_getptset_2d(V(s,:),dist,a);
        succ=tz_validatetredge(V,E,s,pts(2,:));    
        iter=iter+1;
    end
    
    if iter<=maxiter
        degrees(s-1)=deg;
        degrees(end+1)=1;
        V(end+1,:)=pts(2,:);
        size(V,1)
        E(end+1,:)=[s,size(V,1),dist];
    else
        V=orgV;
        E=orgE;
        degrees=1;
    end
end

tz_plotgraph(V,E,[],'x','-',[])

figure
E=tz_mstree(V,'eu',1);
E=E';