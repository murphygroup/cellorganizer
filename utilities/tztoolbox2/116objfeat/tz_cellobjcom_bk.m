function cellobjcom=tz_cellobjcom_bk(classes,cellidcs,post,ncluster,featset,...
    combobj,flindex,cdindex,combcoords,combobjects)
%TZ_CELLOBJCOM_BK Unknown.

%function cellobjcom=tz_cellobjcom(classes,cellidcs,post,ncluster,featset,combobj,flindex,cdindex,combcoords,combobjects)
%OVERVIEW: 
%   Convert object-level features into cell-level features
%PARAMETERS:
%   classes - class indices for n objects nx1
%   cellidcs - cell indices for n objects nx1
%   post - cluster indices for n objects nx1 or nxk
%   ncluster - cluster number, usually it equals to cluster number for post
%   featset - the set of features for return
%       '#' - feature seperator
%       '$' or '!'- cluster marker
%       'num' 'featssum' 'featsmean' 'frac' 'mst' 
%   combobj - p-D features for n objects nxp
%   flindex - fluorescence fraction index, 1x1(just index) or nx1(already features)
%   cdindex - cof dist index,  1x1(just index) or nx1(already features)
%   combcoords - combined cofs for objects
%   combobjects - combined objects
%RETURN:
%   cellobjcom - 1xnclass cell array. 
%       The size for each element depends on featset:
%           'objnum','fluofrac': ncellxk
%           'all': ncellx(2k+k*p) or ncellx(2k+k*(p-1))
%DESCRIPTION:
%   multiclass cell features(MCF): {cells of calss1} {cells of class2} ... {cells of classn}
%
%HISTORY:
%   03-MAR-2004 Initial Write TINGZ
%   19-APR-2004 Modified TINGZ
%       - add tz_findclass
%       - add featset
%   11-MAY-2004 Modified TINGZ
%       - remove cell labels
%   17-MAY-2004 Modified TINGZ
%       - add comments
%       - add mean features
%   26-MAY-2004 Modified TINGZ
%       - add ncluster==-1
%   02-JUN-2004 Modified TINGZ
%       - add cdindex
%   06-JUN-2004 Modified TINGZ
%       - feature combination
%   02-JUL-2004 Modified TINGZ
%       - add combojects

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
caclass=tz_findclass(classes);

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
                case 'fluofrac'
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
             
                if featsets{fi}(1)=='$' |featsets{fi}(1)=='!'
                    
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
                        case 'fluofrac'
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


