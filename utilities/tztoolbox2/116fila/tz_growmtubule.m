function [mts,bidx]=tz_growmtubule(length,n,p,mu,sigma)
%TZ_GROWMTUBULE Under construction.

%function tz_growmtubule(length,n)

theta=2*pi/n;

for i=0:n-1
    mts{i+1}=[(0:length)',zeros(length+1,1)];
%     bidx{i+1}(1)=1;
%     k=2;
%     for j=2:length
%         bend=binornd(1,p,1,1);
%         if bend==1
%             
%             bidx{i+1}(k)=j;
%             k=k+1;
%             prevpt=mts{i+1}(j-1,:);
%             curpt=mts{i+1}(j:end,:);
%             ba=normrnd(mu,sigma,1,1);
%             while ~(ba<pi/2 & ba>-pi/2)
%                 ba=normrnd(mu*(bidx{i+1}(k)-bidx{i+1}(k-1)),sigma,1,1);
%             end
%             curpt=tz_rotate_2d([curpt(:,1)-prevpt(1),curpt(:,2)-prevpt(2)],ba);
%             mts{i+1}(j:end,:)=[curpt(:,1)+prevpt(1),curpt(:,2)+prevpt(2)];
%         end
%         
%     end
%     bidx{i+1}(k)=length+1;
%     mts{i+1}=[mts{i+1}(:,1),spline(mts{i+1}(bidx{i+1},1),mts{i+1}(bidx{i+1},2),mts{i+1}(:,1))];
    mts{i+1}=tz_rotate_2d(mts{i+1},theta*i);
end