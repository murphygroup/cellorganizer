function [features,names]=tz_cellobjfeats(combclass,combcellidx, ...
    combobj,nclass,ncells,objfeatnames)
%TZ_CELLOBJFEATS Calclulate cell-level features from objects.
%   FEATURES = TZ_CELLOBJFEATS(COMBCLASS,COMBCELLIDX,COMBOBJ,
%       NCLASS,NCELLS,OBJFEATNAMES)
%   
%   [FEATURES,NAMES] = TZ_CELLOBJFEATS(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [features,names]=tz_cellobjfeats(combclass,combcellidx,combobj,nclass,ncells,objfeatnames)
%OVERVIEW:
%   calculate cell-level features from objects
%PARAMETERS:
%   combclass - class index
%   combcellidx - cell index
%   nclass - number of classes
%       -1 - automaticly find it in combclass
%   ncells - number of cells in each class
%   combobj - object-level features
%   objfeatnames - object-level feature names
%RETURN:
%   features - cell-level features MCF
%   names - cell-level feature names
%DESCRIPTION:
%
%HISTORY:
%   27-MAY-2004 Initial write TINGZ
%   25-Mar-2004 Modified TINGZ
%       - Change sum to var

if nclass==-1
    caclass=tz_findclass(combclass);
    nclass=length(caclass);
end

if isempty(ncells)
    ncells=-ones(1,nclass);
end

%nclass=max(combclass);
nobjfeats=size(combobj,2);

for i=1:nclass
    if ncells(i)==-1
        cacell=tz_findclass(combcellidx(combclass==i));
        ncell=length(cacell);
    else
        ncell=ncells(i);
        for k=1:ncell
            cacell{k}=k;
        end
    end
    
    features{i}=zeros(ncell,nobjfeats*2+1);
    for j=1:ncell
        objfeats=combobj(combcellidx==cacell{j}(1) & combclass==i,:);
        if ~isempty(objfeats)
            features{i}(j,1:nobjfeats)=mean(objfeats,1);
            if size(objfeats,1)>1
                features{i}(j,nobjfeats+1:end-1)=var(objfeats);
            else
                features{i}(j,nobjfeats+1:end-1)=zeros(1,nobjfeats);
            end
            features{i}(j,end)=size(objfeats,1);
        end
    end
    
end

if isempty(objfeatnames)
    names={};
else
    for i=1:length(objfeatnames)
        names{i}=['average_' objfeatnames{i}];
    end
    
    for i=1:length(objfeatnames)
        names{nobjfeats+i}=['variance_' objfeatnames{i}];
    end
    
    names{end+1}='objnum';
end
