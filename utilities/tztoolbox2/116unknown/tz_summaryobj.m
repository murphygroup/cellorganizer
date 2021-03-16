function [summaryobj,cellobj]=tz_summaryobj(objfeat)

%function [summaryobj,cellobj]=tz_summaryobj(objfeat)
%
%OVERVIEW:
%   
cellidx=objfeat(:,end);
cellnum=max(cellidx);

for i=1:cellnum
    cellobj{i}=objfeat(find(cellidx==i),1:(end-1));
    objnum(i)=size(cellobj{i},1);
end

wwpvalue=zeros(cellnum,cellnum)+NaN;

% for i=1:cellnum
%     for j=(i+1):cellnum
%         if (objnum(i)>100 & objnum(j)>100)
%             wwpvalue(j,i)=tz_wwtest2(cellobj{i},cellobj{j});
%         end
%     end
%     i
% end

summaryobj=struct('objnum',objnum,'wwpvalue',wwpvalue);



