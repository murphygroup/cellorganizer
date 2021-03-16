function E=tz_mstree_l(X)
%TZ_MSTREE Minimal spanning tree for large size data.
%   TREE = TZ_MSTREE_L(X) is similar with TZ_MSTREE. But it only works for
%   Euclidean distance. And it is much faster than TZ_MSTREE on data with
%   large size. TREE and X are the same as those for TZ_MSTREE.
%   
%   See also TZ_MSTREE

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nsample=size(X,1);

pset=2:nsample;
pnewset=1;
pupdate=pnewset;
E=[];
oldindex=[];
olddist=[];
oldupdate=[];
forwoard=1; 
alter=0;
while(~isempty(pset))
    if length(pset)==1
        test=1;
    end
    %if(length(pset)>=length(pnewset))
    length(pset);
    if length(pset)>1
        [index, distance] = nn_search(X(pset,:), ...
            nn_prepare(X(pset,:)), X(pupdate,:), 1);
        index=pset(index);
        
        allindex=[oldindex,index];
        alldist=[olddist;distance];
        allupdate=[oldupdate,pupdate];
        
        [mindist,minj]=min(alldist);  
        mini=allindex(minj);
               
        E=[E;[allupdate(minj),mini,mindist]];
        pnewset=[pnewset,mini];
        pset(pset==mini)=[];
        
        oldindex=allindex(allindex~=mini);
        olddist=alldist(allindex~=mini);
        pupdate=[allupdate(allindex==mini),mini];
        oldupdate=allupdate(allindex~=mini);
    else
        [index, distance] = nn_search(X(pnewset,:), ...
            nn_prepare(X(pnewset,:)), X(pset,:), 1);
        E=[E;pset,pnewset(index),distance];
        pset=[];
    end
%     else
%         if forwoard=1
%             forword=0;
%             alter=1;
%         else
%             alter=0;
%         end
%         [index, distance] = nn_search(X(pnewset,:), nn_prepare(X(pnewset,:)), X(pupdate,:), 1);
%         index=pnewset(index);
%         
%         allindex=[oldindex,index];
%         alldist=[olddist,distance];
%         allupdate=[oldupdate,pupdate];
%         
%         [mindist,minj]=min(alldist);  
%         mini=allindex(minj);
%                
%         E=[E;[allupdate(minj),mini,mindist]];
%         pnewset=[pnewset,mini];
%         pset(pset==mini)=[];
%         
%         oldindex=allindex(allindex~=mini);
%         olddist=alldist(allindex~=mini);
%         pupdate=[pupdate(allindex==mini),mini];
%         oldupdate=pupdate(allindex~=mini);
%end
    
%     if alter==1
%         pupdate=pset;
%         oldindex=[];
%         olddist=[];
%         oldupdate=[];
%     end
end

E=E';
