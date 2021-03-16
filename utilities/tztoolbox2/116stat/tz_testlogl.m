function logl=tz_testlogl(objcom,objnum,testobj,mixp,testnum)

%function logl=tz_testlogl(objcom,objnum,testobj,mixp,testnum)
%
%OVERVIEW:
%   calculate the log-likelihood test data from the multinomial model
%PARAMETERS:
%   objcom - training data, mxn Matrix for n object compositions for m cells
%   objnum - training data, pdf of object number
%   testobj - objects for calculating likelihood
%   mixp - for calculating likelihood of 2-mixture models
%   testnum - object numbers of testing set. 
%             testnum - 0, it will be calculated from test features
%                       -1, no object numbers
%RETURN:
%   logl - log likelihood, a matrix
%DESCRIPTION:
%   this function can also calculate the log-likelihood of 2 mixture
%
%HISTORY:
%   14-MAR-2004 Initial write TINGZ
%   02-SEP-2004 Modified TINGZ
%       - add comments

if ~exist('mixp','var')
    mixp={};
end

if isempty(objnum)
    testnum==-1;
end

if testnum==-1
    objnum=[];
end

if ~exist('testnum','var') | testnum==0
    testnum=sum(testobj,2);
end

nobj=size(testobj,1);
nclass=size(objcom,1);

for i=1:nobj
    for j=1:nclass
        p=tz_normobjcom(objcom(j,:)+tz_normobjcom(testobj(i,:)));
        if ~isempty(objnum)
            logl(i,j)=tz_mnomlogl(testobj(i,:),p,tz_discrdens(objnum{j},50,testnum(i)));
        else
            logl(i,j)=tz_mnomlogl(testobj(i,:),p);
        end
    end
    
    for j=1:length(mixp)
        mix=mixp{j}(1,:);
        trs=mixp{j}(2,:);
        p=tz_normobjcom(objcom(mix(1),:)+objcom(mix(2),:)*trs(1)/trs(2)+tz_normobjcom(testobj(i,:)));
        logl(i,j+nclass)=tz_mnomlogl(testobj(i,:),p);
        if ~isempty(objnum)
            p1=tz_normobjcom(objcom(mix(1),:)+tz_normobjcom(testobj(i,:)));
            p2=tz_normobjcom(objcom(mix(2),:)+tz_normobjcom(testobj(i,:)));
            alpha=tz_fitmnom(testobj(i,:),p1,p2);
            x1=round(alpha*sum(testobj(i,:)));
            x2=sum(testobj(i,:))-x1;
            logl(i,j+nclass)=logl(i,j+nclass)+log(tz_discrdens(objnum{mix(1)},500,x1))+...
                log(tz_discrdens(objnum{mix(2)},500,x2));
        end
    end
end

