function [cellfeat,cellnums] = tz_objfeat2cellfeat(objfeats,featset,weights)
%TZ_OBJFEAT2CELLFEAT Convert object features to cell feature.
%   CELLFEAT = TZ_OBJFEAT2CELLFEAT(OBJFEATS,FEATSET,WEIGHTS)
%   
%   [CELLFEAT,CELLNUMS] = TZ_OBJFEAT2CELLFEAT(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function cellfeat = tz_objfeat2cellfeat(objfeats,featset)
%OVERVIEW
%   convert object features to cell feature
%PARAMETERS
%   objfeats - cell array
%   featset - feature set: 'avg' 'var' 'num'
%RETURN
%   cellfeat - vector
%   cellnums - number of cells in each class
%DESCRIPTION
%   
%HISTORY
%   15-May-2005 Initial write TINGZ
%SEE ALSO
%   tz_cellobjfeats

if ~exist('weights','var')
    weights=[];
end

combobjfeats=[];
objnums=[];
ws=[];
cellnums=zeros(1,length(objfeats));
for i=1:length(objfeats)
    for j=1:length(objfeats{i})
        combobjfeats=[combobjfeats;objfeats{i}{j}];
        objnums=[objnums;size(objfeats{i}{j},1)];
        if ~isempty(weights)
            ws=[ws;weights{i}{j}];
        end
        if ~isempty(objfeats{i}{j})
            cellnums(i)=cellnums(i)+1;
        end
    end
    
end

cellfeat=[];
for i=1:length(featset)
    switch featset{i}
    case 'avg'
        cellfeat=[cellfeat,mean(combobjfeats,1)];
    case 'var'
        if size(combobjfeats,1)>1
            cellfeat=[cellfeat,var(combobjfeats)];
        else
            cellfeat=[cellfeat,zeros(size(combobjfeats))];
        end
    case 'num'
        cellfeat=[cellfeat,sum(objnums)];
    case 'wvg'  %weigthed average
        cellfeat=[cellfeat,ws'*combobjfeats/sum(ws)];
    case 'wvr'
        for j=1:size(combobjfeats,2)
            cellfeat=[cellfeat,ml_wmoment(combobjfeats(:,j)',ws',2)];
        end
    case 'wtm'
        for j=1:size(combobjfeats,2)
            cellfeat=[cellfeat,ml_wmoment(combobjfeats(:,j)',ws',3)];
        end
    end
end
