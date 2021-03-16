function [objcom,objnum]=tz_trainobjcom(cellobjcom,sel)
%TZ_TRAINOBJCOM Summarize object composition.
%   OBJCOM = TZ_TRAINOBJCOM(CELLOBJCOM,SEL) returns a vector of the fraction 
%   of objects  of each type. CELLOBJCOM is a cell array of object number
%   matrices. Each row of the matrix CELLOBJCOM{I} is the numbers of 
%   different types of objects in a cell of class I. SEL is a vector of
%   indices for cell selection. If an element of SEL is negative, the
%   corresponding cell is excluded.
%   
%   [OBJCOM,OBJNUM] = TZ_TRAINOBJCOM(...) also returns the distribution of
%   total object number in a cell. OBJNUM has two rows, the first row 
%   contains the numbers of objects in a cell and the second row contains
%   the numbers of cells.

%   14-MAR-2004 Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University



nclass=length(cellobjcom);
sel2=sel;
sel(sel<0)=[];
for i=1:nclass
    cellobjcom{i}(-sel2(sel2<0),:)=[];
    
    if ~isempty(sel)
        objcom(i,:)=sum(cellobjcom{i}(sel,:),1);
        objnum{i}=[sum(cellobjcom{i}(sel,:),2),ones(size(cellobjcom{i}(sel,:),1),1)];
    else
        objcom(i,:)=sum(cellobjcom{i},1);
        objnum{i}=[sum(cellobjcom{i},2),ones(size(cellobjcom{i},1),1)];
    end
end
