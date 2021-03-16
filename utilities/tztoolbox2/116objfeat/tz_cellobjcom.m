function cellobjcom=tz_cellobjcom(classes,cellidcs,post,ncluster, ...
    featset,combobj,flindex,cdindex,combcoords,combobjects)
%TZ_CELLOBJCOM Convert object-level features into cell-level features.
%   CELLOBJCOM = TZ_CELLOBJCOM(CLASSES,CELLIDCS,POST,NCLUSTER,FEATSET,
%       COMBOBJ,FLINDEX,CDINDEX,COMBCOORDS,COMBOBJECTS) returns a cell
%   array of cell-level features, which are from the clustering labels
%   POST (1-N coding or just labels). CLASSES is a vector of class 
%   labels and CELLIDS is a vector of cell indices. NCLUSTER is the number 
%   of clusters. It is usually equal to cluster number implied by POST.
%   FEATSET is a string of feature sets. Each feature set follows after
%   '#'. See TZ_CALCOBJBASEDSET for the details of the available feature 
%   sets. If there is a '$' or '!' (old design) between '#' and feature
%   set name, the features will be calculated on each cluster.
%   FLINDEX specifies the object fluorescence fractions. If it is a
%   vector, FLINDEX is the vector of object fluorescence fractions. 
%   If it is a scalar, FLINDEX is the column index of fluorescence
%   fraction in COMBOBJ. CDINDEX specifies the vector of object COF
%   to cell COF distance. It can also be a vector or scalar and
%   similar with FLINDEX. COMBCOORS is a matrix of object center
%   coordinates. COMBOBJECTS is a one-level cell array of all objects.
%   The returned value CELLOBJCOM has MCF structure. COMBOBJ, FLIINDX, 
%   CDINDEX, COMBCOORDS and COMBOBJECTS could be empty if they are not
%   needed for the specified feaure set. The size of each feature matrix 
%   depends on feature sets.
%
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   03-MAR-2004 Initial Write T. Zhao
%   19-APR-2004 Modified T. Zhao
%       - add tz_findclass
%       - add featset
%   11-MAY-2004 Modified T. Zhao
%       - remove cell labels
%   17-MAY-2004 Modified T. Zhao
%       - add comments
%       - add mean features
%   26-MAY-2004 Modified T. Zhao
%       - add ncluster==-1
%   02-JUN-2004 Modified T. Zhao
%       - add cdindex
%   06-JUN-2004 Modified T. Zhao
%       - feature combination
%   02-JUL-2004 Modified T. Zhao
%       - add combojects
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 10
    error('Exactly 10 arguments are required')
end

%string to cell array
featsets=tz_parsefeatset(featset);

%object types
if size(post,2)>1
    [m,objidcs]=max(post');
    objidcs=objidcs';
else
    objidcs=post;
end

%nx3 [object classes, object cell indices, object types]
objs=[classes,cellidcs,objidcs];

%find class. caclass is a cell array of [class label,sample indices]
caclass=ml_findclass(classes);

%number of classes
nclass=length(caclass);

%if number of clusters is -1, it will be set to the maximum of cluster labels
if ncluster==-1
   ncluster=max(objidcs);
end

% if ~exist('featset','var')
%     featset='#objnum';
% end

% if ~exist('combobj','var')
%     combobj=[];
% end

%index of fluorescence fraction.
%If the given parameter is a scalar, 
%it is supposed to be the index of object features
if exist('flindex','var')
    if(length(flindex)==1)
        tmpindex=flindex;
        flindex=combobj(:,flindex);
        combobj(:,tmpindex)=[];
    end
end

if exist('cofdistindex','var')
    if(length(cdindex)==1)
        tmpindex=cdindex;
        cdindex=combobj(:,cdindex);
        combobj(:,tmpindex)=[];
    end
end

%Use cell indices to calculate object numbers if combobj is empty
if isempty(combobj)
    combobj = cellidcs;
end

cellobjcom={};

for i=1:nclass
    classi=objs(objs(:,1)==i,:);
    ncell=max(classi(:,2));
    objcom=[];
    for j=1:ncell
        cellj=classi(classi(:,2)==j,:);
        
        if ~isempty(cellj)
            allobjs=combobj(objs(:,1)==i & objs(:,2)==j,:);
            cellobj=[];
            for fi=1:length(featsets)
                switch featsets{fi}
                case {'featsum','featmean'}
                    if isempty(allobjs)
                        inputobj=zeros(1,size(combobj,2));
                    else
                        inputobj=allobjs;
                    end
                case 'objnum'
                    inputobj=allobjs;
                case {'fluofrac','totalfluo','objsize'}
                    inputobj=flindex(objs(:,1)==i & objs(:,2)==j);
                case {'mst','mst2'}
                    inputobj=combcoords(objs(:,1)==i & objs(:,2)==j,:);
                case 'fmst'
                    inputobj=[allobjs,combcoords(objs(:,1)==i & objs(:,2)==j,:)];
                case 'readymst'
                    load /home/tingz/matlab/tz_objfunction/readyfeats/readymst.mat
                    [classi(1,1),cellj(1,2)]
                    inputobj=mstfeats{classi(1,1)}(cellj(1,2),:);
                case {'pixelnum','fluosum'}
                    inputobj=combobjects(objs(:,1)==i & objs(:,2)==j);
                otherwise
                    inputobj=[];
                end
                
                %%%%%%%%for debugging only%%%%%%%%%
                [i j]
%                 dbg=tz_getenv('dbg');
%                 if dbg{1}==1
%                     if i==9 & j==31
%                         'stop'
%                     end
%                 end
                %%%%%%%%%%%%%%%%%%%%%%
                
                cellobj=[cellobj,tz_calcobjbasedset(inputobj,featsets{fi})];
             
                if featsets{fi}(1)=='$' |featsets{fi}(1)=='!' ...
                        | featsets{fi}(1)==':'
                    
                    for k=1:ncluster
                        allobjsclst=combobj(objs(:,1)==i & objs(:,2)==j & objs(:,3)==k,:);
                        switch featsets{fi}(2:end)
                        case {'featsum','featmean'}
                            if isempty(allobjsclst)
                                inputobj=zeros(1,size(combobj,2));
                            else
                                inputobj=allobjsclst;
                            end
                        case 'objnum'
                            inputobj=allobjsclst;
                        case {'fluofrac','totalfluo','objsize'}
                            inputobj=flindex(objs(:,1)==i & objs(:,2)==j & objs(:,3)==k);
                        case {'mst','mst2'}
                            inputobj=combcoords(objs(:,1)==i & objs(:,2)==j & objs(:,3)==k,:);
                        case 'fmst'
                            inputobj=[allobjsclst,combcoords(objs(:,1)==i & objs(:,2)==j & objs(:,3)==k,:)];
                        case {'pixelnum','fluosum'}
                            inputobj=combobjects(objs(:,1)==i & objs(:,2)==j & objs(:,3)==k);
                        otherwise
                            inputobj=[]; 
                        end
                        cellobj=[cellobj,tz_calcobjbasedset(inputobj,featsets{fi}(2:end))];
                    end
                end
            end
            
            objcom=[objcom;cellobj];
        end
    end
    cellobjcom{i}=objcom;
end


