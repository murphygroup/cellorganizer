function [ncvcm,avgacc]=tz_summarymix2(b)

%function tz_summarymix2(b)
%
%OVERVIEW:
%   Summarize 2 mixture results
%PARAMETERS:
%   b - 2mixture cross validataion results
%RETURN
%   ncvcm - average confusion matrix
%   avgacc - average accuracy
%DESCRIPTION:
%   
%HISTORY:


nfold=length(b);
npair=size(b{1},1);
ncvcm=zeros(npair,npair);

for i=1:nfold
    for j=1:npair
        for k=1:size(b{i},2)
            index=b{i}(j,k);    
            if index~=0
                ncvcm(j,index)=ncvcm(j,index)+1;
            end
        end
    end
end

avgacc=sum(diag(ncvcm))/sum(ncvcm(:));
