function ncm = tz_summarytestclassif_pslid(y,pred_y,tr_y)
%TZ_SUMMARYTESTCLASSIF_PSLID Obsolete.

%function ncm = tz_summarytestclassif_pslid(y,pred_y,tr_y)
%OVERVIEW
%   
%PARAMETERS
%   y - 
%   pred_y - 
%   try - 
%RETURN
%   ncm - 
%DESCRIPTION
%   
%HISTORY
%   24-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_summarytestclassif_pslid','ml_summarytestclassif'));

ntestclass=max(y);
ntrainclass=max(tr_y);

for i=1:ntestclass
    for j=1:ntrainclass
        index=find(y==i);
        if isempty(index)
            ncm(i,j)=0;
        else
            ncm(i,j)=sum(double(pred_y(index)==j));
        end
    end
end
