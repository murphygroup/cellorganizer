function n=tz_countsame(ind1,ind2)
%TZ_COUNTSAME Count the number of values that occur in both vectors.
%   N = TZ_COUNTSAME(IND1,IND2) returns the number of values that occur
%   in both vectro IND1 and IND2.

%   27-NOV-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

m=max(max(ind1),max(ind2));
hist1=hist(ind1,m);
hist2=hist(ind2,m);
n=sum(hist1&hist2);

%??? n=length(intersect(ind1,ind2)) ???

% n=0;
% while ~isempty(ind1) & ~isempty(ind2)
%     
%     if ind1(1)==ind2(1)
%         n=n+1;
%         ind1(1)=[];
%         ind2(1)=[];
%         continue
%     end
%     
%     if ind1(1)>ind2(1)  %swap
%         tmp=ind1;
%         ind1=ind2;
%         ind2=tmp;
%     end
%     
%     dind=ind1-ind2(1);
%     ind1(dind<0)=[];
%     dind(dind<0)=[];
%     ind2(1)=[];
%     if ~isempty(dind)
%         if dind(1)==0
%             n=n+1;
%         end
%     end
% end