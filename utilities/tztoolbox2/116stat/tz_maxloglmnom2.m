function logl=tz_maxloglmnom2(objcom,objnum,testobj,mixp)
%TZ_MAXLOGLMNOM2 Unknown.
%   LOGL = TZ_MAXLOGLMNOM2(OBJCOM,OBJNUM,TESTOBJ)
%   
%   LOGL = TZ_MAXLOGLMNOM2(OBJCOM,OBJNUM,TESTOBJ,MIXP)
%   
%   See also

%   19-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function logl=tz_maxloglmnom2(objcom,objnum,testobj,mixp)
%
%OVERVIEW:
%   calculate maximum log likelihood for each mixp
%PARAMETERS:
%   objcom - training object composition
%   objnum - object number distribution
%   testobj - testing objects 
%   mixp - which 2 mixtures from objcom
%RETURN:
%   logl - loglikehood matirx
%DESCRIPTION:
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments

if ~exist('mixp','var')
    mixp={};
end

ncell=size(testobj,1);
nclass=size(objcom,1);

for i=1:ncell
    for j=1:length(mixp)
        mix=mixp{j}(1,:);
        
        p1=tz_normobjcom(objcom(mix(1),:)+tz_normobjcom(testobj(i,:)));
        p2=tz_normobjcom(objcom(mix(2),:)+tz_normobjcom(testobj(i,:)));
        alpha=tz_fitmnom2_ga(testobj(i,:),p1,p2);
        if(alpha==0 | alpha==1)
            logl(i,j)=-Inf;
        else
            x1=round(alpha*sum(testobj(i,:)));
            x2=sum(testobj(i,:))-x1;
            logl(i,j)=tz_mnomlogl(testobj(i,:),alpha*p1+(1-alpha)*p2)+log(tz_discrdens(objnum{mix(1)},500,x1))+...
                log(tz_discrdens(objnum{mix(2)},500,x2));
        end
        
    end
end

