function namelists = tz_classiflist(pred_y,y,setnames)
%TZ_CLASSIFLIST Obsolete

%function namelists = tz_classiflist(pred_y,y,setnames)
%OVERVIEW
%   
%PARAMETERS
%   pred_y - predicted label
%   y - true label, must be 1...n
%   setnames - name of classes
%RETURN
%   namelists - cell array of name lists of each true class
%DESCRIPTION
%   
%HISTORY
%   25-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_classiflist','ml_classiflist'));

nclass=max(y);

for i=1:nclass
    pred_yi=pred_y(y==i);
    namelists{i}=setnames(pred_yi);
end
