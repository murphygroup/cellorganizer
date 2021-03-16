function [synimg,psimg]=tz_logisticimg(celledge,ra,nucedge,para,N)
%TZ_LOGISTICIMG Synthesize an image by logistic model.
%   SYNIMG = TZ_LOGISTICIMG(CELLEDGE,RA,NUCEDGE,PARA,N)
%   
%   [SYNIMG,PSIMG] = TZ_LOGISTICIMG(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [synimg,psimg]=tz_logisticimg(celledge,ra,nucedge,para,N)
%
%OVERVIEW:
%   generate an image by logistic model
%PARAMETERS:
%   celledge - cell edge image
%   ra - orientation of the cell
%   nucedge - nucleus edge image
%   para - paramerets of the logistic model
%   N - max intensity
%RETURN:
%   synimg - synthesized image
%   psimg - probability image
%DESCRIPTION:
%   
%HISTORY:
%   ??-OCT-2004 Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add comments

synimg=zeros(size(celledge));
psimg=synimg;


celldist=bwdist(celledge);
nucdist=bwdist(nucedge);

cellimg=imfill(celledge,'hole');
nucimg=imfill(nucedge,'hole');

codemask=double(cellimg)-double(nucimg)-double(celledge);

[I,J]=find(codemask==1);
ind=sub2ind(size(codemask),I,J);

if isempty(ra)    
    dists=[nucdist(ind),celldist(ind)];
else
    polcodes=tz_cellpolcode(celledge,ra,nucedge,{},0);
    dists=polcodes(:,1:2);
end

ps=tz_evallogistic([ones(size(dists,1),1),dists],para);
psimg(ind)=ps;

for i=1:length(ind)
    synimg(ind(i))=binornd(N,ps(i),1,1);
end

